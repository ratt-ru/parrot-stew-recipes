#!/usr/bin/env python
import os, re
import os.path
import shutil
import pickle
import traceback
import fnmatch
import datetime
import warnings
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import numpy.ma as ma
from numpy import fft
from scipy.ndimage import gaussian_filter, fourier_gaussian
from scipy.interpolate import interp1d
from astropy import units as u
from astropy.coordinates import SkyCoord, SpectralCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time, TimeDelta
from astropy.coordinates import FK5
from astropy.units.quantity import Quantity
from astropy.visualization import time_support
from astropy.table import QTable
from casacore.tables import table

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from typing import Optional, Any, List
from regions import Regions
import pandas as pd
import click

# boilerplate for supporting click wrappers around python cabs
from scabha.schema_utils import clickify_parameters
from stimela.kitchen.cab import Cab
from omegaconf import OmegaConf
schemas = OmegaConf.load("parrot-cabs.yml")
@click.group()
def cli():
    pass


def get_freq_axis(hdr):
    for faxis in range(1, hdr['NAXIS']+1):
        if hdr[f"CTYPE{faxis}"] == "FREQ":
            return faxis
    else:
        return None

def stack_xarray_cube(images: List[str], cube: str, ms: str, verbose: bool = False, cadence=1, include_freq_axis: bool = False):
    



def stack_time_cube(images: List[str], cube: str, ms: str, verbose: bool = False, cadence=1, include_freq_axis: bool = False):
    print(f"will process {len(images)} snapshot images")
    hdr = fits.open(images[0])[0].header
    # get MS timestamps
    ms_tms = sorted(set(table(ms).getcol("TIME")))
    t0 = ms_tms[0]
    dt = ms_tms[1] - t0
    grid = np.arange(t0, ms_tms[-1] + dt/2, dt)
    print(f"MS has {len(ms_tms)} unique timestamps. Total grid of {len(grid)} timestamps with interval {dt}s")
    # get freq axis and copy to axis 4
    axiskw = 'NAXIS', 'CTYPE', 'CUNIT', 'CRPIX', 'CRVAL', 'CDELT'
    faxis = get_freq_axis(hdr)
    if faxis is not None:
        if include_freq_axis:
            if faxis != 4:
                for kw in axiskw:
                    hdr[f"{kw}4"] = hdr[f"{kw}{faxis}"]
            hdr['NAXIS'] = 4
        else:
            cfreq = SpectralCoord(hdr[f"CRVAL{faxis}"], hdr[f"CUNIT{faxis}"])
            dfreq = SpectralCoord(hdr[f"CDELT{faxis}"], hdr[f"CUNIT{faxis}"])
            hdr['FREQ0MHZ'] = (cfreq - dfreq/2).to_value(u.MHz)
            hdr['FREQ1MHZ'] = (cfreq + dfreq/2).to_value(u.MHz)
            faxis = None
    if faxis is None:
        # no frequency info, drop axis 4
        for kw0 in axiskw:
            kw = f"{kw0}4"
            if kw in hdr:
                del hdr[kw]
        hdr['NAXIS'] = 3
        # add keywords
    # form up time axis
    hdr['NAXIS3'] = len(grid)
    hdr['CTYPE3'] = 'TIME'
    hdr['CUNIT3'] = 'sec'
    hdr['CRPIX3'] = 1
    hdr['CRVAL3'] = 0
    hdr['CDELT3'] = dt
    # 
    cubedata = np.zeros((len(grid), hdr['NAXIS2'], hdr['NAXIS1']), np.float32)
    max_delta = 0
    for nimg, image in enumerate(images):
        ff = fits.open(image)
        tm = Time(ff[0].header['DATE-OBS'])
        # convert to MS time representation
        t = tm.mjd * 24 * 3600 - t0
        it = round(t/dt)
        offset = t - it*dt
        max_delta = max(max_delta, abs(offset))
        if verbose:
            print(f"timeslot {it} for {tm} (image {nimg}): offset {offset}s")
        if abs(offset) > dt/2:
            raise ValueError(f"timeslot {it} for {tm} (image {nimg}) is too far from a time grid point ({offset}s)")
        if it < 0 or it >= len(grid):
            raise ValueError(f"timeslot {it} for {tm} (image {nimg}) out of range")
        # paste in
        cubedata[it, :, :] = ff[0].data[0, 0, :, :]
        # deallocate
        ff.close()
    # reshape for freq axis
    if faxis is not None:
        cubedata = cubedata.reshape([1]+list(cubedata.shape))

    print(f"max time offset was {max_delta}s")
    fits.writeto(cube, cubedata, hdr, overwrite=True)
    print(f"wrote time cube to {cube}")


def convolve_image(image, outimage, size_arcsec=None, size_pix=None, size_sec=0, use_fft=False, ncpu: int=0):
    if image == outimage:
        return

    time0 = datetime.datetime.now()
    def logtime():
        return str(datetime.datetime.now() - time0)

    ff = fits.open(image, memmap=True)
    hdr = ff[0].header

    if size_pix is None:
        if size_arcsec is None:
            size_arcsec = 3600*(hdr['BMAJ'] + hdr['BMIN'])/2
        
        scale = 3600*hdr['CDELT1']
        size_pix = abs(size_arcsec / scale) / 2.355 # convert FWHM to sigma

    size_timeslots = size_sec / hdr['CDELT3']

    sigma = (size_timeslots, size_pix, size_pix)

    slicer = tuple([0]*(hdr['NAXIS'] - 3) + [slice(None)]*3)

    print(f"{logtime()} reading image")        
    # make copy of image
    img = ff[0].data[slicer]
    # make mask of time axis
    print(f"{logtime()} calculating mask and time weights")        
    mask = (img==0).all(axis=(1,2))
    # make weights vector by convoving the inverse of the mask
    weight = np.ones_like(mask, float)
    weight[mask] = 0
    weight = gaussian_filter(weight, sigma=sigma[0], mode='constant')

    if use_fft and not size_pix:
        output = img
        print(f"{logtime()} padding and transposing array")
        padding = int(sigma[0]*3)
        pimg = np.zeros(img.shape[1:] + (2*padding + img.shape[0],), img.dtype)
        pimg[:, :, padding:-padding] = img.transpose((1, 2, 0))
        # reshape array in memory
        print(f"{logtime()} computing rfft")
        imgfft = fft.rfft(pimg, axis=2)
        print(f"convolving with kernel {sigma}")
        fourier_gaussian(imgfft, sigma[0], axis=2, output=imgfft)
        print(f"computing rifft")
        imgfft = fft.irfft(imgfft, axis=2)
        output[...] = imgfft[...,padding:-padding].transpose((2,0,1))
    else:
        print(f"{logtime()} convolving {image} with kernel {sigma}")
        # convolve image to output
        output = img
        gaussian_filter(img.copy(), sigma, output=output, mode='constant')

    print(f"{logtime()} applying time weights")
    output /= weight[..., np.newaxis, np.newaxis]
    output[mask, :, :] = 0

    print(f"{logtime()} writing {outimage}")
    ff.writeto(outimage, overwrite=True)


def extract_fits_metadata(images, timestamps_file, beams_file):
    print(f"extracting timestamps and beam info from {len(images)} images")
    timestamps = []
    beams = []
    for i, image in enumerate(images):
        ff = fits.open(image)
        tm = Time(ff[0].header['DATE-OBS'])
        timestamps.append(tm.mjd)
        beams.append(list(ff[0].header[key] for key in ('BMIN', 'BMAJ', 'BPA')))
    pickle.dump(timestamps, open(timestamps_file, "wb"))
    pickle.dump(beams, open(beams_file, "wb"))

cat_tab = None
fluxes = None

def _export_regions(numbers: List[int], outfile: str, default_color: str = "yellow"):
        with open(outfile, "wt") as rf:
            rf.write("""# Region file format: DS9 CARTA 2.0.0
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 
fk5\n""")
            for num in numbers:
                label = cat_tab['labels'][num]
                lbl = f"{num} ({label})" if label else str(num)
                color = cat_tab['color'][num] or default_color
                pos = cat_tab['pos'][num]
                rad = cat_tab['R'][num]
                flux = fluxes[num]
                rf.write(f"""circle({pos.ra.deg}, {pos.dec.deg}, 20") # color={color} width=2 sep={rad} flux={flux} text={{{lbl}}}\n""")
        print(f"wrote {outfile}")


wcs = None
img = None
img_xy = {}
freq0 = None
freq1 = None
chunk_dir = None
stdbox_r = None
filelabel = ""
std_flag_threshold = None
sigma_threshold = 4
subtract_curves = None
output_dir = None
beam_info = None

timestamps = None
highlight_timestamps = []
plot_title_format = None

def _extract_curve(isrc: int, srcnum: int):

    what = "light curve" if not subtract_curves else "baseline-subtracted lightcurve"   
    src = cat_tab[srcnum]
    flux = fluxes[srcnum]
    x, y = img_xy[srcnum]
    radius = src['R'].to_value(u.deg)
    print(f"extracting {what} #{isrc} for source ID {srcnum}, flux {flux*1e+6:.2f} uJy at {x}, {y} (dist {radius:.2f}deg)")

    try:
        coord = src['pos']
        lc_file = f"{chunk_dir}/src{srcnum:06d}{filelabel}"

        lightcurve = img[x,y,:]
        stddev = img[x-stdbox_r:x+stdbox_r+1,y-stdbox_r:y+stdbox_r+1,:].std(axis=(0,1))
        # mask gaps
        mask = lightcurve==0

        # mask timeslots where stddev exceeds median
        medstd = np.median(stddev[~mask])
        if std_flag_threshold:
            stdmask = stddev > std_flag_threshold*medstd
            mask |= stdmask
        else:
            stdmask = None
            
        if mask.all():
            print(f"lightcurve {srcnum} is empty, is there any data?")
            return srcnum, None

        lightcurve = np.ma.masked_array(lightcurve, mask)
        stddev = np.ma.masked_array(stddev, mask)

        lightcurve_raw = lightcurve
        # get beam gain if available, and correct the lightcurve for it
        if beam_info:
            bgain = beam_info.get_time_variable_beamgain(coord, mjds=timestamps)[0]
            bgain1 = beam_info.get_time_variable_beamgain(coord, mjds=timestamps, spi=2)[0]
            bgain2 = beam_info.get_time_variable_beamgain(coord, mjds=timestamps, spi=4)[0]
            # lightcurve = lightcurve - flux * (bgain - bgain.mean()) / bgain.mean() 
            lightcurve = lightcurve_raw - flux * (bgain / bgain.mean() - 1)
            lightcurve1 = lightcurve_raw - flux * (bgain1 / bgain1.mean() - 1)
            lightcurve2 = lightcurve_raw - flux * (bgain2 / bgain2.mean() - 1)
            # print(f"{srcnum} {(bgain / bgain.mean() - 1)}")
        else:
            bgain = None

        # subtract if asked
        if subtract_curves:
            sub_lc_file = lc_file.replace(output_dir, subtract_curves)
            sub_lc, _, _ = pickle.load(open(f"{sub_lc_file}.p", "rb"))
            lightcurve -= sub_lc

        if sigma_threshold:
            significant = np.where(abs(lightcurve) > sigma_threshold*stddev)
        else:
            significant = []

        peak = lightcurve.max()
        peak_std = stddev[lightcurve.argmax()]
        maxvar = abs(lightcurve).max() / flux
        stdvar = lightcurve.std() 
        relstdvar = stdvar / flux

        # make lightcurve plot
        ra_s = coord.ra.to_string(u.hour, precision=1)
        dec_s = coord.dec.to_string(u.degree, precision=1)
        ra_s0 = coord.ra.to_string(u.hour, precision=0)
        dec_s0 = coord.dec.to_string(u.degree, precision=0)

        nts = len(lightcurve)
        fig = plt.figure(figsize=(20,4))
        sz = 1 * 3000./nts
        sz = max(sz, min(sz, 1), 0.5)
        timeslot_axis = np.arange(nts) 
        # plot base lightcurves
        plt.errorbar(timeslot_axis, lightcurve*1e+6, yerr=stddev*1e+6, 
                        fmt='.', mec='navy', mfc='navy', ecolor='deepskyblue', ms=sz, elinewidth=sz, zorder=10)
        # plot beam-corrected lightcurves, if given
        if bgain is not None:
            plt.errorbar(timeslot_axis, lightcurve1*1e+6, yerr=stddev*1e+6, 
                            fmt='.', mec='navy', mfc='navy', ecolor='deepskyblue', ms=sz, elinewidth=sz, zorder=10)
            plt.errorbar(timeslot_axis, lightcurve2*1e+6, yerr=stddev*1e+6, 
                            fmt='.', mec='navy', mfc='navy', ecolor='deepskyblue', ms=sz, elinewidth=sz, zorder=10)
        # plot signficant excursions
        if significant:
            plt.errorbar(timeslot_axis[significant], lightcurve[significant]*1e+6, yerr=stddev[significant]*1e+6, 
                            fmt='.', mec='red', mfc='red', ecolor='pink', ms=sz, elinewidth=sz, zorder=10)
        # plot masked events
        if stdmask is not None:
            plt.plot(timeslot_axis[stdmask], np.zeros_like(timeslot_axis)[stdmask], 'x',
                                mec='purple', mfc='purple', zorder=5)
        plt.xlabel("timeslot number")
        ylabel = "Flux offset, uJy"
        if freq0 is not None:
            ylabel += f" ({freq0:.0f}-{freq1:.0f} MHz)"
        if subtract_curves:
            ylabel += " (baseline subtracted)"
        plt.ylabel(ylabel)
        plt.xlim(0, nts-1)
        plt.hlines([0], 0, nts-1, colors='gray', ls='dotted', zorder=2, lw=1)

        if bgain is not None:
            plt.plot(lightcurve_raw*1e+6, marker='.', mec='green', mfc='green', ms=sz, zorder=8)
            plt2 = plt.twinx()
            plt2.plot(bgain, '.', ms=sz, zorder=9, color='grey')
            plt2.plot(bgain1, '.', ms=sz,  zorder=9, color='grey')
            plt2.plot(bgain2, '.', ms=sz,  zorder=9, color='grey')
            plt2.set_ylabel("Power beam gain")

        ylim = plt.ylim()
        plt.vlines(np.arange(0, nts, 20), *ylim, colors='lightgray', ls='dotted', zorder=1, lw=1)
        plt.vlines(np.arange(0, nts, 100), *ylim, colors='gray', ls='dotted', zorder=2, lw=1)

        with warnings.catch_warnings(record=True) as w:
            titleprops = dict(
                label=f"{srcnum} ({src['labels']})" if src['labels'] else f"Source {srcnum}",
                ra_str=ra_s, dec_str=dec_s,
                flux=flux, flux_ujy=flux*1e+6,
                medstd=medstd, medstd_ujy=medstd*1e+6,
                peak=peak, peak_ujy=peak*1e+6,
                peak_std=peak_std, peak_std_ujy=peak_std*1e+6,
                radius=radius,
                maxvar=maxvar,
                relstdvar=relstdvar,
                vm=stdvar/medstd
            )
            if len(w):
                print(f"source {srcnum} generated {len(w)} warnings")

        title = plot_title_format or "Source {label}, {ra_str} {dec_str}, {flux_ujy} uJy, r={radius:.2f}deg, " + \
                    "median std={medstd_ujy:.2f}uJy, variation max {maxvar:.2%}, rms {relstdvar:.2%}, vm {vm:.2f}"
        plt.title(title.format(**titleprops), color=src['color'] or 'black')

        # highlight interesting events in lightcurve of source 0
        if highlight_timestamps and srcnum == 0:
            for ts, color, label in highlight_timestamps:
                plt.vlines([ts], ylim[0], lightcurve[ts]*1e+6, colors=color or 'green', lw=1)

        lc_plot_file = f"{chunk_dir}/lc{srcnum:06d}-{ra_s0}{dec_s0}-r{radius:.2f}d{filelabel}.png"

        # disable UserWarning from AutoDateLocator
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")        
            # add time axis
            with time_support(simplify=True):
                ax = plt.twiny()
                ax.set_xlabel('Time (UTC)')
                t0, t1 = timestamps[0].to_datetime(), timestamps[-1].to_datetime()
                # print(t0, t1)
                ax.set_xlim(t0, t1)
                # locator = mdates.MinuteLocator(interval=15)
                locator = mdates.AutoDateLocator(minticks=24)
                formatter = mdates.ConciseDateFormatter(locator)
                ax.xaxis.set_major_locator(locator)
                ax.xaxis.set_major_formatter(formatter)
            fig.set_tight_layout(True)

            plt.savefig(lc_plot_file, dpi=300, bbox_inches='tight', pad_inches=0)
            plt.close()

        # make power spectrum plot
        plt.figure(figsize=(20,2))
        f = fft.fft(lightcurve)
        f[0] = 0
        f = fft.fftshift(f)
        plt.plot(abs(f), '.', ms=sz)
        plt.xlabel("Fourier mode")
        plt.ylabel("Abs power")
        plt.title(title)
        ps_plot_file = f"{chunk_dir}/powerspectrum{srcnum:06d}{filelabel}.png"
        plt.savefig(ps_plot_file, dpi=300)
        plt.close()

        # make stats dict
        dfdict = {}
        # convert to kosher types (for QTable(meta) to take it)
        for key, value in zip(src.keys(), src.values()):
            if isinstance(value, Quantity):
                if key.startswith("e"):
                    value = value.to_value(u.arcsec)
                elif key.startswith("flux"):
                    value = value.to_value(u.Jy)
                else:
                    continue
            elif isinstance(value, np.str_):
               value = str(value) 
            elif isinstance(value, SkyCoord) or isinstance(value, np.ndarray):
                continue
            dfdict[key] = value

        dfdict.update(x=x, y=y, flux=flux, 
            ra=coord.ra.deg, dec=coord.dec.deg,
            maxvar=maxvar,
            stdvar=stdvar,
            medstd=medstd,
            meanstd=stddev.mean(),
            lc_plot=lc_plot_file, ps_plot=ps_plot_file, lc_file=lc_file)
        table_meta = dfdict.copy()

        if freq0 is not None:
            table_meta.update(freq0=Quantity(freq0, u.MHz), freq1=Quantity(freq1, u.MHz))
        table_meta['pos'] = coord

        # write lightcurve to file
        cols = OrderedDict(time=timestamps, flux=lightcurve, flux_err=stddev)
        units = OrderedDict(flux="Jy", flux_err="Jy")
        if bgain is not None:
            cols.update(beamgain=bgain, flux_raw=lightcurve_raw)
            units.update(flux_raw="Jy")
        tab = QTable(list(cols.values()), names=list(cols.keys()), units=units, meta=table_meta)

        try:
            tab.write(f"{lc_file}.ecsv", overwrite=True)
        except Exception as exc:
            with open(f"{lc_file}.err", "wt") as ef:
                ef.write(traceback.format_exc())
                ef.write("\n")
                ef.write(f"{timestamps}\n")
                ef.write(f"{lightcurve}\n")
                ef.write(f"meta: {table_meta}\n")
            pickle.dump((lightcurve, timestamps, stddev), open(f"{lc_file}.p"))

        return srcnum, pd.DataFrame(dfdict, index=[srcnum])

    except Exception as exc:
        traceback.print_exc()
        print(f"extracting light curve {srcnum}: failed with exception, see above")
        return srcnum, None

sources = None

def extract_light_curves(cube, catalog, outdir, statsfile=None, regfile=None, nsrc=100, chunk_size=100, 
                         # use this field as the flux
                         fluxcols=None,
                         # matches interesting timestamps from text file
                         interesting_timestamps=None,
                         select_labels=['*'],    # include matching labels. '*' includes unlabelled sources
                         # additional criteria for unlabeled sources 
                         within=None,             # ...within this distance of center (Quantity string)
                         srctype=None,            # ...with this type (P/G)
                         minflux=None,               # ...above this flux (Quantity string)
                         # plot title, otherwise default title is used
                         plot_title=None, #
                         # other options
                         beaminfo=None,           # BeamInfo file for power beam, to enable beam correction
                         subtract=None,           # baseline lightcurves to subtract
                         sigma=4,                 # threshold for marking lightcurve excusrsions
                         flag_excess_std=None,    # if set, time bins with rms> this*median_rms are flagged
                         output_file_label="",
                         stdbox=50,               # size of box for rms estimates 
                         default_color="yellow",  # default colour for unlabelled sources
                         ncpu=4):
    # make these global, so that multiprocessing can pickle them
    global wcs, img, chunk_dir, stdbox_r, filelabel, std_flag_threshold, sigma_threshold, sources, \
            subtract_curves, output_dir, beam_info, timestamps, highlight_timestamps, plot_title_format
    std_flag_threshold = flag_excess_std
    stdbox_r = stdbox
    filelabel = output_file_label
    sigma_threshold = sigma
    subtract_curves = subtract
    output_dir = outdir
    plot_title_format = plot_title

    if beaminfo:
        print(f"Loading power beam info from {beaminfo}")
        beam_info = pickle.load(open(beaminfo, "rb")) 

    # load catalog
    global cat_tab
    cat_tab = QTable.read(catalog)

    for fluxcol in fluxcols:
        if fluxcol not in cat_tab.colnames:
            raise RuntimeError(f"table {catalog} foes not contain column '{fluxcol}'")
    global fluxes
    # make Ncolumns,Nrows array of flux columns
    fluxes = np.array([cat_tab[fluxcol].to_value(u.Jy) for fluxcol in fluxcols])
    # get first non-zero flux along each row
    fluxes = fluxes[np.argmax(fluxes!=0, axis=0), range(fluxes.shape[1])]
    # get field centre
    centre = cat_tab.meta['centre']
    # get labels
    labels = [set(lbl.split(":")) for lbl in cat_tab['labels']]

    # read image cube
    ff = fits.open(cube)
    hdr = ff[0].header
    global freq0, freq1
    freq0 = hdr.get("FREQ0MHZ", None)
    freq1 = hdr.get("FREQ1MHZ", None)
    # get wcs and reduce to 2 axes
    wcs = WCS(hdr)
    while len(wcs.array_shape) > 2:
        wcs = wcs.dropaxis(len(wcs.array_shape) - 1)
    img = ff[0].data.T

    # build array of MJD timestamps
    t0 = Time(hdr['DATE-OBS'], format='fits') + TimeDelta(hdr['CRVAL3'], format=hdr['CUNIT3'])
    dt = TimeDelta(hdr['CDELT3'], format=hdr['CUNIT3'])
    timestamps = t0 + dt*np.arange(hdr['NAXIS3'])

    if interesting_timestamps:
        highlight_timestamps = []
        for line in open(interesting_timestamps, "rt").readlines():
            if '#' in line:
                line = line[:line.index("#")]
            line = line.strip()
            if line:
                match = re.match("(?P<time>[^\s]+)(\s+(?P<color>[^\s]+))?(\s+(?P<label>[^\s]+))?", line)
                if not match:
                    raise RuntimeError(f"can't parse entry in {interesting_timestamps}: {line}")
                # get matching timestamp, if any
                if match['time'].startswith("20"):
                    mjd = Time(match['time'], format='isot')
                else:
                    mjd = Time(match['time'], format='mjd')
                delta_t = abs((mjd - timestamps).sec)
                imin = np.argmin(delta_t)
                dt_sec = delta_t[imin]
                if dt_sec < 16:
                    print(f"event at {match['time']} matched within {dt_sec:.02f}s")
                    highlight_timestamps.append((imin, match['color'], match['label']))
                else:
                    print(f"event at {match['time']} not matched to this lightcurve span ({dt_sec:.02f}s)")

    # build list of sources to plot
    all_sources = list(range(len(cat_tab)))
    sel_sources = set()

    for sel in select_labels:
        # select by matcing label
        if sel != "*":
            subset = [num for num, lbls in enumerate(labels) if any(fnmatch.fnmatch(lbl, sel) for lbl in lbls)]
            sel_sources.update(subset)
            print(f"label '{sel}' selects {len(subset)} sources")
        else:
            sources = all_sources
            if minflux:
                print(f"selecting sources brighter than {minflux}")
                minflux = Quantity(minflux).to_value("Jy") 
                sources = [num for num in sources if fluxes[num] >= minflux]
                print(f"selection by flux leaves {len(sources)} sources")
            if within:
                print(f"selecting sources within {within} of centre")
                within = Quantity(within).to_value("deg")
                radius = cat_tab['R'].to_value("deg")
                sources = [num for num in sources if radius[num] <= within]
                print(f"selected {len(sources)} sources within {within:.2f} degrees")
            if srctype:
                print(f"selecting sources of type {srctype}")
                types = cat_tab['type']
                sources = [num for num in sources if fnmatch.fnmatch(types[num], srctype)]
                print(f"selection by type leaves {len(sources)} sources")
            sel_sources.update(sources)

    sources = sorted(sel_sources)
    print(f"{len(sources)} selected for lightcurves")
    if not sources:
        return

    # compute pixel positions and eliminate same-pixel and out-of-bounds sources
    img_xp, img_yp = wcs.world_to_pixel(cat_tab["pos"][sources])
    img_xp = np.round(img_xp).astype(int)
    img_yp = np.round(img_yp).astype(int)
    # make global dict of isrc -> (x, y)
    global img_xy
    img_xy = {isrc: (x, y) for isrc, x, y in zip(sources, img_xp, img_yp)}

    # go through list of positions in reverse, and form dict from {(x,y):src}, so that for multiple sources at the same (x,y), 
    # the lowest-numbered source ends up in the dictionary. Also eliminate out-of-bounds pixels
    posdict = {(x, y): isrc for isrc, x, y in list(zip(sources, img_xp, img_yp))[::-1]
                if x>=0 and y>=0 and x<img.shape[0] and y<img.shape[1]}

    sources = list(posdict.values())[::-1]
    print(f"eliminating non-unique and out-of-bounds pixel positions leaves {len(sources)} sources")

    # print the ones eliminated
    print(f"  (eliminated: {sorted(set(sel_sources) - set(sources))})")

    if nsrc is not None and len(sources) > nsrc:
        sources = sources[:nsrc]
        print(f"  restricting to first {len(sources)}")

    nsrc = len(sources)

    # now start extracting
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.mkdir(outdir)

    # dump region file
    if regfile is not None:
        _export_regions(sources, regfile, default_color=default_color)

    stats_df = None

    for ichunk, i0 in enumerate(range(0, nsrc, chunk_size)):
        i1 = min(i0 + chunk_size, nsrc) 

        # make output directory
        chunk_dir = f"{outdir}/src-chunk{ichunk:05d}"
        if not os.path.exists(chunk_dir):
            os.mkdir(chunk_dir)

        _export_regions(sources[i0:i1], 
                        f"{chunk_dir}/regions-chunk{ichunk:05d}.reg",
                        default_color=default_color)

        if ncpu > 1:
            with ProcessPoolExecutor(ncpu) as pool:
                # submit each iterant to pool
                futures = [pool.submit(_extract_curve, i0+i, isrc) for i, isrc in enumerate(sources[i0:i1])]
                        
                results = [f.result() for f in futures]
        else:
            results = []
            for i, isrc in enumerate(sources[i0:i1]):
                results.append(_extract_curve(i0+i, isrc))

        stats = [st for _, st in results if st is not None]

        if stats:
            if stats_df is None:
                stats_df = pd.concat(stats, ignore_index=True)
            else:
                stats_df = pd.concat([stats_df]+stats, ignore_index=True)

    if statsfile:
        pickle.dump(stats_df, open(statsfile, "wb"))
    

model_hdus = model_center_freq = model_delta_freq = model_freq0 = model_freq1 = None
model_wcs = None
model_shape = None
model_slice = None
freqgrids = None
freqdeltas = None
finest_freqgrid = None

def _extract_single_spectrum(lcfile: str):
    # no, give spectrum per model bin
    lctab = QTable.read(lcfile)
    freq0 = lctab.meta['freq0'].to_value(u.Hz)
    freq1 = lctab.meta['freq1'].to_value(u.Hz)
    # get pixel position
    x, y = model_wcs.world_to_pixel(lctab.meta['pos'])
    x, y = round(float(x)), round(float(y))
    # outside of model image? bug out
    if x < 0 or y < 0 or x >= model_shape[0] or y >= model_shape[1]:
        return None
    # extract overlap with finest model grid
    if freq0 > finest_freqgrid[-1] or freq1 < finest_freqgrid[0]:
        raise RuntimeError("lightcurve frequencies don't overlap model frequencies")
    # [i0:i1] is the frequency grid within the interval of the lightcurve then 
    i0 = np.where(finest_freqgrid >= freq0)[0][0]
    i1 = np.where(finest_freqgrid <= freq1)[0][-1] + 1
    print(f"lightcurve {lcfile} frequency {freq0*1e-6:.2f} to {freq1*1e-6:.2f} MHz corresponds to model grid [{i0}:{i1}]")
    #
    freqs = finest_freqgrid[i0:i1]
    model = np.zeros_like(freqs)
    # now scan all relevant models
    for hdus, grid, delta in zip(model_hdus, freqgrids, freqdeltas):
        wh0 =  np.where(grid + delta/2 >= freqs[0])[0]
        wh1 =  np.where(grid - delta/2 <= freqs[-1])[0]
        if not len(wh0) or not len(wh1):
            continue
        # j0:j1 is the overlap interval of the modelset with the frequency grid
        j0, j1 = wh0[0], wh1[-1]+1
        # read off model flux
        modelflux = [hdu.data[(model_slice)+(y,x)] for hdu in hdus[j0:j1]]
        # interpolate onto grid
        print(f"using model grid [{j0}:{j1}], flux is {modelflux}")
        model += interp1d(grid[j0:j1], modelflux, "nearest", fill_value="extrapolate")(freqs)
    print(f"Estimated model flux {model.mean()}, source finder flux was {lctab.meta['flux']}")
    print(f"Model spectrum is {model}")
    return model.mean()


@cli.command("ems")
@clickify_parameters(Cab(**schemas.cabs['extract-model-spectrum']))
def extract_model_spectrum(lightcurves: List[str], modelsets: List[List[str]], ncpu: int=1):
    global model_hdus
    model_hdus = [[fits.open(filename, memmap=True)[0] for filename in models] for models in modelsets]
    # get header of first model, assume the grid is uniform throughout
    hdr0 = model_hdus[0][0].header
    ifreq = get_freq_axis(hdr0)
    if ifreq is None:
        return RuntimeError("no FREQ axis in first model")
    # get WCS from first model
    global model_wcs
    model_wcs = WCS(hdr0)
    while len(model_wcs.array_shape) > 2:
        model_wcs = model_wcs.dropaxis(len(model_wcs.array_shape) - 1)
    global model_shape, model_slice
    model_shape = (hdr0['NAXIS1'], hdr0['NAXIS2'])
    model_slice = (0,)*(hdr0['NAXIS'] - 2)
    # get frequency span of each modelset, and figure out which is the finest grid,
    # as that's the one we'll interpolate the models onto
    global freqgrids, freqdeltas, finest_freqgrid
    freqgrids = []
    freqdeltas = []
    finest_freqgrid = []
    for hdus in model_hdus:
        freqs = np.array([hdu.header[f"CRVAL{ifreq}"] for hdu in hdus])
        if len(freqs) > len(finest_freqgrid):
            finest_freqgrid = freqs
        freqgrids.append(freqs)
        freqdeltas.append(np.array([hdu.header[f"CDELT{ifreq}"] for hdu in hdus]))
    if ncpu > 1:
        with ProcessPoolExecutor(ncpu) as pool:
            # submit each iterant to pool
            futures = [pool.submit(_extract_single_spectrum, lc) for lc in lightcurves]
            results = [f.result() for f in futures]
    else:
        resutls = [_extract_single_spectrum(lc) for lc in lightcurves]

if __name__ == '__main__':
    cli()




    