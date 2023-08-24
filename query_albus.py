#!/usr/bin/env python
# A basic python script that tracks a specified position on the 
# sky over the time range from START_TIME to END_TIME from
# a specific location on the Earth's surface.

# The output is a text file giving Slant Tec (STEC) and
# ionosphere rotation measure (RM) as a function of time


import os
import time
import sys
import matplotlib
import subprocess
matplotlib.use("agg") # display should not be set

import math

def __run_rinex_predict(RA, DEC, LONG, LAT, HEIGHT, OBJECT, DATA_DIR, MAX_DIST, START_TIME, END_TIME):
    import MS_Iono_functions as iono 
    
    os.system('date')
    process_start = time.time()
    startime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("getGPSIONO Start at %s" % startime)

    # Moon experiment
    RED_TYPE = 'RI_G01'
    TIME_STEP = 300

    NUM_PROCESSORS = 1
    # Note: we set NUM_PROCESSORS = 1 as we are getting data from Geosciences 
    # Australia, which seems to have difficulty responding to a number of ftp
    # requests being received in parallel 
    # After the GPS data have been collected the system will increase the
    # number of processors for the final ionosphere modelling
    iono. process_ionosphere(time_step=TIME_STEP,
                             object=OBJECT,
                             Ra=RA,
                             Dec=DEC,
                             Lat=LAT,
                             Long=LONG,
                             Height=HEIGHT,
                             start_time=START_TIME,
                             end_time=END_TIME,
                             max_dist=MAX_DIST,
                             processing_option=RED_TYPE,
                             do_serial=0,
                             num_processors=NUM_PROCESSORS,
                             gps_data_directory=DATA_DIR)

    os.system('date')
    endtime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("getGPSIONO End at %s" % endtime)
    print (' ')
    process_end = time.time()
    duration = (process_end - process_start)/60.0
    print("getGPSIONO Total run time: %7.2f minutes" % duration)
    print("current directory is", os.getcwd())

def __checkmake_outputdir():
    testoutdir = os.environ.get("ALBUS_TESTCASE_OUTPUT", "/output_dir")
    if os.path.exists(testoutdir):
        outputhiddendir = os.path.join(testoutdir, 'datadumps')  
        if not os.path.exists(outputhiddendir):
            os.mkdir(outputhiddendir)
    else:
        raise RuntimeError(f"'{testoutdir}' must be a valid directory -- are you running inside Docker, check mount path?")
    return outputhiddendir

def run_albus_predict(ms: str,  output_dir: str, location: str = "MeerKAT", field: int = 0):
    from casacore.tables import table
    from astropy.coordinates import SkyCoord, FK5, EarthLocation
    from astropy.time import Time
    import astropy.units as u

    # get field info
    field_tab = table(f"{ms}::FIELD")
    OBJECT = field_tab.getcol("NAME", field)[0]
    ra, dec = field_tab.getcol("PHASE_DIR", field)[0,0]
    field_tab.close()
    pos = SkyCoord(ra, dec, unit=u.rad, frame=FK5)
    RA = pos.ra.to_string(u.hour, sep=":", precision=2)
    DEC = pos.dec.to_string(u.deg, sep=":", precision=2)

    # get earth location
    if location: 
        loc = EarthLocation.of_site(location)
    else:
        raise RuntimeError("must read location from MS, but this is not implemented")

    LON = loc.lon.to_string(u.deg, sep=":", precision=2)
    LAT = loc.lat.to_string(u.deg, sep=":", precision=2)
    HEIGHT = str(loc.height.to_value(u.m))

    # get time
    timecol = table(ms).getcol("TIME")
    START_TIME = Time(timecol.min()/(24*3600), format='mjd').to_datetime().strftime("%Y/%m/%d %H:%M:%S")
    END_TIME = Time(timecol.max()/(24*3600), format='mjd').to_datetime().strftime("%Y/%m/%d %H:%M:%S")

    print(f"Direction is {OBJECT} at {RA} {DEC}")
    print(f"Location is {LON} {LAT} {HEIGHT}")
    print(f"Time range is {START_TIME} {END_TIME}")

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if not os.path.exists(f"{output_dir}/tecs"):
        os.mkdir(f"{output_dir}/tecs")

    
    args = ["docker", "run", 
        "-v", "/home/oms/projects/ALBUS_ionosphere/gfzrnx_2.0-8219_lx64:/optsoft/bin/gfzrnx",
        "-v", "/home/oms/projects/ALBUS_ionosphere:/albus_waterhole",
        "-v", f"{os.path.abspath(__file__)}:/query_albus.py",
        "-v", f"{os.path.abspath(output_dir)}:/output_dir",
        "--workdir", "/output_dir", 
        "--rm", "--user", f"{os.getuid()}:{os.getgid()}",
        "albus:latest",
        "/query_albus.py",
        RA, DEC, LON, LAT, HEIGHT, START_TIME, END_TIME, "target", "/output_dir/rinex"]
    
    print("running:", " ".join(args))
    subprocess.check_call(args)

if __name__ == "__main__":
    print(sys.argv)

    RA, DEC, LONG, LAT, HEIGHT, START_TIME, END_TIME, OBJECT, DATA_DIR = sys.argv[1:]

    MAX_DIST = 350E3
    __run_rinex_predict(RA, DEC, LONG, LAT, float(HEIGHT), OBJECT, DATA_DIR, MAX_DIST, START_TIME, END_TIME)
