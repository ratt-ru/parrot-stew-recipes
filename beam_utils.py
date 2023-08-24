import pickle
from typing import Any, List, Optional
from dataclasses import dataclass
import os.path
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u
from scipy.ndimage import spline_filter, map_coordinates

@dataclass
class PowerBeam(object):
    """Power beam info"""
    I: np.ndarray      # Stokes I beam, of shape NFREQ x NDEG x NDEG
    deg: np.ndarray    # coordinates in beam
    freq: np.ndarray   # frequencies

def mdv_beams_to_power_beam(mdv_beams: str, power_beam: str):
    """
    Converts MdV's npz beamset into a Stokes I power beam
    """
    mdv = np.load(mdv_beams)
    bm = mdv['beam']
    degs = mdv['margin_deg']
    freqs = mdv['freq_MHz']*1e+6
    i0 = len(degs) // 2
    xx2 = abs(bm[0,-1,:,:,:])**2
    yy2 = abs(bm[3,-1,:,:,:])**2
    # normalize by beam gain at centre, since the bandpass will do that for us
    xx2norm = xx2 / xx2[:,i0,i0][:,np.newaxis,np.newaxis]
    yy2norm = yy2 / yy2[:,i0,i0][:,np.newaxis,np.newaxis]
    inorm = (xx2norm + yy2norm)/2
    hdu = fits.PrimaryHDU(inorm[:,::-1,:])
    hdr = hdu.header
    hdr['CRPIX1'] = i0+1
    hdr['CRPIX2'] = i0+1
    hdr['CRPIX3'] = 1
    hdr['CRVAL1'] = 0
    hdr['CRVAL2'] = 0
    hdr['CRVAL3'] = freqs[0]
    hdr['CDELT1'] = degs[1] - degs[0]
    hdr['CDELT2'] = degs[1] - degs[0]
    hdr['CDELT3'] = freqs[1] - freqs[0]
    hdr['CTYPE1'] = 'X'
    hdr['CTYPE2'] = 'Y'
    hdr['CTYPE3'] = 'FREQ'
    hdr['CUNIT1'] = 'deg'
    hdr['CUNIT2'] = 'deg'
    hdr['CUNIT3'] = 'Hz'
    hdu.writeto(f"{power_beam}", overwrite=True)
    print(f"saved Stokes I beam to {power_beam}")

def derive_power_beam(cube: str, images: List[str], beaminfo: str, 
                      outcube: Optional[str] = None, nband: Optional[int] = None, 
                      power_beam: str ='beam/MeerKAT_U_band_StokesI.npz'):
    """
    Derives a power beam across the specified number of subbands 
    (same as the bands in the FITS cube, if nband is not supplied),
    and assigns weights to each subband using the WSCVWSUM value from 
    the image list.
    If outcube is given, saves beam-corrected cube as well.
    """
    # load beam info
    pbf = fits.open(power_beam)   # axes are FREQ,Y,X
    I = pbf[0].data.transpose((0,2,1))  # now FREQ,X,Y
    hdr = pbf[0].header
    x0, y0 = hdr['CRPIX1']-1, hdr['CRPIX2']-1
    nx, ny = hdr['NAXIS1'], hdr['NAXIS2']
    delta = hdr['CDELT1']
    # beam frequencies
    pbfreq_step = hdr['CDELT3']
    pbfreq = hdr['CRVAL3'] + (np.arange(hdr['NAXIS3']) - hdr['CRPIX3'] + 1)*pbfreq_step
    minfreq, maxfreq = pbfreq[0] - pbfreq_step/2, pbfreq[-1] + pbfreq_step/2 

    # cube
    cubef = fits.open(cube)
    cubehdr = cubef[0].header 
    ra0 = cubehdr['CRVAL1']
    dec0 = cubehdr['CRVAL2']
    # find FREQ axis
    for faxis in range(1, cubehdr["NAXIS"]+1):
        if cubehdr[f"CTYPE{faxis}"] == "FREQ":
            print(f"Cube FREQ axis is {faxis} (1-based)")
            break
    else:
        raise RuntimeError("Input cube does not have a FREQ axis")
    # numpy freq axis
    faxis_np = cubehdr["NAXIS"] - faxis
    # work out freq grid
    cube_nfreq = cubehdr[f'NAXIS{faxis}'] 
    cube_dfreq = cubehdr[f'CDELT{faxis}']
    cube_freqs = cubehdr[f'CRVAL{faxis}'] + \
        (np.arange(cube_nfreq) - cubehdr[f'CRPIX{faxis}'] + 1)*cube_dfreq
    print(f"Cube frequencies are {cube_freqs}")
    # slicing object to select freq planes
    cube_slice = [0]*(cubehdr['NAXIS']-2) + [None, None]
    if len(images) != cube_nfreq:
        raise RuntimeError(f"cube has {cube_nfreq} frequencies, but {len(images)} images given")    

    if not nband:
        nband = cube_nfreq
    # subband frequencies
    freq_step = (maxfreq-minfreq)/nband
    freqs = minfreq + (np.arange(nband) + 0.5)*freq_step

    # populate sum-weights from FITS images
    sumweights = np.zeros_like(pbfreq)
    image_band_masks =[]
    for freq, image in zip(cube_freqs, images):
        # freq range of this subimage
        freq0 = freq - cube_dfreq/2
        freq1 = freq + cube_dfreq/2
        # assign sw accordingly
        mask = (pbfreq>=freq0)&(pbfreq<freq1)
        image_band_masks.append(mask)
        sumweights[mask] = fits.open(image)[0].header['WSCVWSUM']

    # compute mean sumweights per band
    band_weights = np.zeros(nband)

    # compute mean beam for each frequency band
    e = np.zeros((nband, nx, ny))
    for band, freq in enumerate(freqs):
        freq0, freq1 = freq - freq_step/2, freq + freq_step/2
        mask = (pbfreq>=freq0)&(pbfreq<freq1)
        band_weights[band] = sumweights[mask].sum()  
        e[band] = (I[mask]*sumweights[mask,np.newaxis,np.newaxis]).sum(axis=0) / band_weights[band]

    # compute mean MFS beam
    emean = (e*band_weights[:,np.newaxis,np.newaxis]).sum(axis=0) / band_weights.sum()

    bi = ObsSpecificBeamInfo(Eband=e, Emean=emean, band_weights=band_weights, freqs=freqs, 
                            x0=x0, y0=y0, delta=delta, ra0=ra0, dec0=dec0)

    pickle.dump(bi, open(beaminfo, "wb"))

    print(f"saved power beam info to {beaminfo}")

    if outcube:
        # make image grid in degrees
        inx, iny = cubehdr['NAXIS1'], cubehdr['NAXIS2']
        ix0, iy0 = cubehdr['CRPIX1'] - 1, cubehdr['CRPIX2'] - 1
        idx, idy = cubehdr['CDELT1'], cubehdr['CDELT2']
        # image grid in degrees
        ixgrid = (np.arange(0, inx) - ix0) * idx
        iygrid = (np.arange(0, iny) - iy0) * idy
        # convert to beam coordinates
        bxgrid = ixgrid / delta + x0
        bygrid = iygrid / delta + y0
        # resample beam to image pixels
        for iband in range(cube_nfreq):
            # make average beam over image band
            Iavg = I[image_band_masks[iband]].mean(axis=0) 
            # resample to image grid
            Iimg = map_coordinates(Iavg, np.meshgrid(bxgrid, bygrid), mode='constant')
            Iimg = Iimg.reshape(inx, iny).transpose()
            # write out
            cube_slice[faxis_np] = iband
            cubef[0].data[tuple(cube_slice)] /= Iimg
        cubef.writeto(outcube, overwrite=True)
        print(f"saved beam-corrected cube as {outcube}")


@dataclass
class ObsSpecificBeamInfo(object):
    """Observation-specific beam info generated by derive_power_beam"""
    Eband: np.ndarray          # per-band power beam
    Emean: np.ndarray          # mean MFS beam
    band_weights: np.ndarray   # per-band weights
    freqs: np.ndarray          # band frequencies
    x0: int       # center pixel of beam
    y0: int
    delta: float  # degrees per pixel
    ra0: float    # field centre in degrees
    dec0: float

    def __post_init__(self):
        self._prefilter_emean = spline_filter(self.Emean)
        self._prefilter_eband = [spline_filter(self.Eband[i,...]) for i in range(len(self.freqs))]

    def get_source_pixel_coodinates(self, srcpos: SkyCoord, mjds: List[float], loc: EarthLocation = None,
                                    signs=(1,1), swap=False):
        """
        Given a sky position and a list of times, derives the in-beam coordinates of the source (in beam pixels)
        """
        centre = SkyCoord(self.ra0*u.deg, self.dec0*u.deg)
        if loc is None:
            loc = EarthLocation.of_site("MeerKAT") 
        # convert positions to alt-az
        frame = AltAz(obstime=Time(mjds, format='mjd'), location=loc)       
        altaz_src = srcpos.transform_to(frame)
        altaz_centre = centre.transform_to(frame)
        # get angle and separation of source w.r.t. centre
        angles = altaz_centre.position_angle(altaz_src)
        seps = altaz_centre.separation(altaz_src)
        # convert to pixel position 
        # confused about angles, but experiments show that 0 is up and +90 is right 
        x = signs[0] * seps.deg * np.sin(angles.rad)
        y = signs[1] * seps.deg * np.cos(angles.rad)
        if swap:
            x, y = y, x
        xp = x / self.delta + self.x0
        yp = y / self.delta + self.y0        
        return np.array([xp, yp])

    def get_time_freq_beamgain(self, srcpos: SkyCoord, mjds: List[float], 
                                loc: EarthLocation = None,
                                signs=(1,1), swap=False):
        xpyp = self.get_source_pixel_coodinates(srcpos, mjds, loc, signs=signs, swap=swap)
        bg_freq = np.zeros((len(self.freqs), len(mjds)))
        for i in range(len(self.freqs)):
            bg_freq[i] = map_coordinates(self._prefilter_eband[i], xpyp, prefilter=True)
        return bg_freq, self.band_weights, xpyp


    def get_time_variable_beamgain(self, srcpos: SkyCoord, mjds: List[float], 
                                        loc: EarthLocation = None, spi: Optional[float] = None):
        xpyp = self.get_source_pixel_coodinates(srcpos, mjds, loc)
        if spi is None:
        # interpolate mean beam
            return map_coordinates(self._prefilter_emean, xpyp, prefilter=True), xpyp
        else:
            bg_freq = np.zeros((len(self.freqs), len(mjds)))
            for i, freq in enumerate(self.freqs):
                bg_freq[i] = map_coordinates(self._prefilter_eband[i], xpyp, prefilter=True)
            spectral_weights = (self.freqs/self.freqs[0])**spi
            weights = self.band_weights*spectral_weights
            return (bg_freq*weights[:,np.newaxis]).sum(axis=0) / weights.sum(), xpyp


