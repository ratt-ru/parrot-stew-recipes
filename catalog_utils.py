from typing import List, Optional, Dict, Tuple, Any
from collections import OrderedDict
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time, TimeDelta
from astropy.coordinates import FK5, concatenate
from astropy.units.quantity import Quantity
from astropy.table import QTable
from astropy.visualization import time_support
from regions import Regions
import numpy as np
import Tigger


def match_catalogs(master_catalog: str, 
                    catalogs: Dict[str, Tuple[str, float, str]], 
                    ra0: float, dec0: float,     
                    max_radius_deg: Optional[float] = None,     # max distance from centre
                    # matches sources from region files
                    interesting_regions: List[str] = [],  # region files for "interesting" sources
                    # matches sources from search box
                    search_box_radec=None,   # if True, separate search in given box of [ra,dec,size,size]
                                            #    (ra, dec are coordinate strings, size is a quantity string)
                    search_box_frame='geocentricmeanecliptic',   # if not None, box is converted to this coordinate frame
                    search_box_label="search box",  # label for matches within this box
                    search_box_minflux=None, #     only consider sources above this flux (Quantity string)
                    search_box_color="green"):
    
    # load model, filter, and convert to coordinates
    num_catalogs = len(catalogs)
    centre = SkyCoord(ra0, dec0, frame=FK5)

    cats = OrderedDict()

    for icat, (label, (catalog, xmatch_arcsec, cat_type)) in enumerate(catalogs.items()):
        model = Tigger.load(catalog)
        coords = SkyCoord(np.array([src.pos.ra for src in model.sources]),
                          np.array([src.pos.dec for src in model.sources]), unit="rad", frame=FK5)
        cats[label] = (model.sources, coords, catalog, xmatch_arcsec, cat_type)  
        print(f"loaded {len(model.sources)} sources from '{cat_type}' {catalog}")
            
    main_icats = [i for i, cat in enumerate(cats.values()) if cat[-1] == "main"]
    print(f"main catalogs are {main_icats}")

    # SPI catalogs
    # cross-match
    for icat, (label, (sources, coords, cat_filename, xmatch_arcsec, cat_type)) in enumerate(cats.items()):
        nsrc = len(sources)
        nums = np.arange(nsrc)
        if icat == 0:
            coords0 = coords
            # table giving source number in each catalog, -1 for not matched
            srcnums = np.full((nsrc, num_catalogs), -1, int)
            # table giving separation w.r.t. each catalog
            srcseps = np.zeros_like(srcnums, float)
            srcnums[:, 0] = nums  # source numbers in first catalog
            # table giving fluxes
            fluxes = np.zeros_like(srcnums, float)
            fluxes[:, 0] = np.array([src.flux.I for src in sources])
            print(f"initializing master catalog from {cat_filename}")
        else:
            # cross-match to master coordinate list
            index, d2d, _ = coords.match_to_catalog_sky(coords0)
            seps = np.array([d.arcsec for d in d2d])
            matched = seps <= xmatch_arcsec
            # update master entries
            srcnums[index[matched], icat] = nums[matched]
            srcseps[index[matched], icat] = seps[matched]
            fluxes[index[matched], icat] = np.array([src.flux.I for src, match in zip(sources, matched) if match])
            print(f"matched {matched.sum()}/{nsrc} sources from {cat_filename} with d={xmatch_arcsec}\"")
            # add rows for missing sources, if catalog is being merged in
            if cat_type == 'main':
                nadd = nsrc - matched.sum()
                # add new objects
                if nadd:
                    print(f"adding {nadd} new sources")
                    newnums = np.full((nadd, num_catalogs), -1, int)
                    newseps = np.zeros_like(newnums, float)
                    newfluxes = np.zeros_like(newnums, float)
                    newnums[:, icat] = nums[~matched]
                    newfluxes[:, icat] = np.array(
                        [src.flux.I for src, match in zip(sources, matched) if not match])
                    coords0 = concatenate((coords0, coords[~matched]))
                    srcnums = np.vstack((srcnums, newnums))
                    srcseps = np.vstack((srcseps, newseps))
                    fluxes  = np.vstack((fluxes, newfluxes))

    nsrc = len(coords0)
    print(f"final catalog contains {nsrc} sources")

    centresep = np.array([coord.separation(centre).deg for coord in coords0])
    if max_radius_deg:
        include = centresep <= max_radius_deg
        coords0 = coords0[include]
        srcnums = srcnums[include]
        srcseps = srcseps[include]
        fluxes = fluxes[include]
        centresep = centresep[include]
        nsrc = len(coords0)
        print(f"selection by R<={max_radius_deg}deg leaves {nsrc} sources")

    # read region files
    interesting_objects = []
    for region_file in interesting_regions:
        regs = Regions.read(region_file, format='ds9')
        for reg in regs:
            label = reg.meta.get('text') or f"OOI{num_match}"
            color = reg.visual.get('edgecolor')
            interesting_objects.append((reg.center, reg.radius, label, color)) 
        print(f"{region_file} contains {len(regs)} objects of interest")

    # add labels and colors from cross-matching interesting regions
    labels = np.full(nsrc, '', dtype=object)
    colors = np.full(nsrc, '', dtype=object)
    UNRANKED = 9999999
    ranks  = np.full(nsrc, UNRANKED)

    # helper function to tag interesting sources
    num_match = 0
    def add_match(idx, label, color):
        nonlocal num_match
        labels[idx] = f"{labels[idx]}:{label}" if labels[idx] else label
        if not colors[idx]:
            colors[idx] = color
        ranks[idx] = min(ranks[idx], num_match)
        num_match += 1

    if interesting_objects:
        print(f"Will try to match {len(interesting_objects)} objects of interest")
        oois = SkyCoord(np.array([pos.ra.deg for pos, _, _, _ in interesting_objects])*u.deg,
                        np.array([pos.dec.deg for pos, _, _, _ in interesting_objects])*u.deg)
        index, d2d, d3d = oois.match_to_catalog_sky(coords0)
        for i, (_, radius, label, color) in enumerate(interesting_objects):
            sep = radius.to_value("arcsec")
            if d2d[i].arcsec <= sep:
                add_match(index[i], label, color)
                print(f"  {i}: matched {label} with source {index[i]}, d={d2d[i].arcsec:.2f}\"")
            else:
                print(f"  {i}: no match for {label} within {sep:.2f}\", nearest source at {d2d[i].arcsec:.2f}\"")

    # search within search box for PSR candidates
    if search_box_radec:
        print(f"Looking for {search_box_label} candidates")
        boxpos = SkyCoord(search_box_radec[0], search_box_radec[1], frame=FK5)
        sep_lon, sep_lat = Quantity(search_box_radec[2]).to_value('arcsec'), Quantity(search_box_radec[3]).to_value('arcsec')
        print(f"Box size is {sep_lon:.2f}\" by {sep_lat:.2f}\"")
        if search_box_frame:
            boxpos = getattr(boxpos, search_box_frame, None)
            if boxpos is None:
                raise RuntimeError(f"invalid search_box_frame '{search_box_frame}'")
            print(f"Box frame is {search_box_frame}")
        if hasattr(boxpos, 'lat'):
            lat0, lon0 = boxpos.lat.arcsec, boxpos.lon.arcsec
        else:
            lat0, lon0 = boxpos.ra.arcsec, boxpos.dec.arcsec
        if search_box_minflux:
            search_box_minflux = Quantity(search_box_minflux).to_value("Jy")
        else:
            search_box_minflux = 0

        candidates = []
        for isrc in range(nsrc):
            if fluxes[isrc, :].max()  > search_box_minflux:
                srcpos = coords0[isrc]
                if search_box_frame:
                    srcpos = getattr(srcpos, search_box_frame)
                if hasattr(srcpos, 'lat'):
                    lat, lon = srcpos.lat.arcsec, srcpos.lon.arcsec
                else:
                    lat, lon = srcpos.ra.arcsec, srcpos.dec.arcsec
                if abs(lat-lat0) < sep_lat and abs(lon-lon0) < sep_lon:
                    print(f"found {search_box_label} candidate at dlat={abs(lat-lat0)/60:.2f}' dlon={abs(lon-lon0)/60:.2f}'")
                    candidates.append((abs(lon-lon0), isrc))
        candidates = [c[1] for c in sorted(candidates)]
        for i, idx in enumerate(candidates):
            add_match(idx, f"{search_box_label}-{i}", search_box_color) 
        print(f"Found {len(candidates)} candidates")

    # now sort by rank, then flux
    sortkey = sorted([(ranks[i], -fluxes[i,main_icats].max(), i) for i in range(nsrc)])
    idx = np.array([key[2] for key in sortkey])
    print(f"sort keys are {' '.join(map(str, idx[:100]))}")

    # organize columns for table
    columns = OrderedDict(
        id        = np.arange(nsrc),
        pos       = coords0[idx],
        R         = centresep[idx],
        cross_id  = srcnums[idx],
        cross_sep = srcseps[idx],
        type      = np.full(nsrc, 'P')
    )
    units = OrderedDict(R=u.deg, cross_sep=u.arcsec)
    NOSPI = 99999.0

    # add flux shape and SPI columns on a per-catalog basis
    for icat, (label, (_, _, _, _, cat_type)) in enumerate(cats.items()):
        columns[f"flux_{label}"] = fluxes[idx, icat]
        columns[f"emaj_{label}"] = np.zeros(nsrc, float)
        columns[f"emin_{label}"] = np.zeros(nsrc, float)
        units[f"flux_{label}"] = u.Jy
        units[f"emaj_{label}"] = u.arcsec
        units[f"emin_{label}"] = u.arcsec
        if cat_type == "spi":
            columns[label] = np.full(nsrc, NOSPI)            

    columns['labels'] = labels[idx]
    columns['color'] = colors[idx]

    # reindex arrays and create table
    tab = QTable(list(columns.values()), 
                    names=list(columns.keys()), 
                    units=units,
                    meta=dict(
                        centre=centre, 
                        catalogs={label: cat[0] for label, cat in catalogs.items()},
                        labels={label: i for i, label in enumerate(catalogs)}
                    ),
    )
    # populate source shape
    nspi = 0
    for icat, (label, (sources, _, cat_filename, _, cat_type)) in enumerate(cats.items()):
        srcnums = tab['cross_id'][:, icat]
        matched = srcnums >= 0
        # source numbers in cross-catalog
        srcnums = srcnums[matched]
        # helper function to extract a nested source attribute, or else return a default value
        def get_attributes(attr1, attr2, default):
            subattr = [getattr(sources[num], attr1, None) for num in srcnums]
            return np.array([getattr(sub, attr2, default) for sub in subattr])
        # match shape info
        emaj = get_attributes('shape', 'ex', 0)
        emin = get_attributes('shape', 'ey', 0) 
        tab[f'emin_{label}'][matched] = Quantity(emin, u.rad)
        tab[f'emaj_{label}'][matched] = Quantity(emaj, u.rad)
        if cat_type == 'spi':
            spi  = get_attributes('spectrum', 'spi', NOSPI) 
            tab[label][matched] = spi
            print(f"{(spi!=NOSPI).sum()} sources got an spi from {cat_filename}")
        # mark source as Gaussian if component is gaussian
        if cat_type == 'main':
            gauss = (emaj>0)|(emin>0)
            tab['type'][matched][gauss] = 'G'
            print(f"{gauss.sum()} sources marked as Gaussian based on {cat_filename}")

    tab.write(master_catalog, overwrite=True)
    print(f"wrote master catalog to {master_catalog}")

def augment_catalog(catalog: str, 
                    augment_catalog: str, 
                    augment_column: str,
                    xmatch_arcsec: float = 6, 
                    coord_column: str = "pos",
                    unmatched_value: Any = '',
                    output_column: str = None,
                    output_catalog: str = None):
    cat0 = QTable.read(catalog)
    cat1 = QTable.read(augment_catalog)
    print(f"supplementary table {augment_catalog} contains {len(cat1)} entries")

    index, d2d, d3d = cat0["pos"].match_to_catalog_sky(cat1[coord_column])
    matched = d2d.arcsec <= xmatch_arcsec
    print(f"cross-matched to {matched.sum()} sources in {catalog} (using d={xmatch_arcsec}\")")

    augcol = cat1[augment_column]
    coldata = [augcol[idx] if match else unmatched_value for match, idx in zip(matched, index)]
    colname = output_column or augment_column
    
    if colname in cat0.colnames:
        cat0[colname] = coldata
    else:
        cat0.add_column(coldata, name=colname)

    cat0.write(output_catalog or catalog, overwrite=True)
    print(f"wrote catalog to {output_catalog or catalog}")






if __name__ == "__main__":
    cols = [[], [], []]
    for line in open("SourceList-MW.txt", "rt").readlines():
        for i, fld in enumerate(line.strip().split()):
            cols[i].append(fld)

    src_id = np.array(cols[0])
    coord = SkyCoord(cols[1], cols[2])

    tab = QTable([src_id, coord], names=("id", "pos"))
    tab.write("SourceList-MW.ecsv", overwrite=True)
