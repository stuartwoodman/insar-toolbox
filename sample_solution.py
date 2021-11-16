import os
import os.path
import subprocess
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

CONF_FILE = r"""
# PyRate configuration file for GAMMA-format interferograms
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Optional ON/OFF switches - ON = 1; OFF = 0

# Coherence masking (PREPIFG)
cohmask:       0

# Orbital error correction (CORRECT)
orbfit:        1

# APS correction using spatio-temporal filter (CORRECT)
apsest:        1

# DEM error (residual topography) correction (CORRECT)
demerror:      1

# Phase Closure correction (CORRECT)
phase_closure:      1

# Optional save of numpy array files for output products (MERGE)
savenpy:       1

# Optional save of incremental time series products (TIMESERIES/MERGE)
savetsincr:    1

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Multi-threading parameters used by correct/stacking/timeseries
# gamma prepifg runs in parallel on a single machine if parallel = 1
# parallel: 1 = parallel, 0 = serial
parallel:      0
# number of processes
processes:     4

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Input/Output file locations
#
# File containing the list of interferograms to use.
ifgfilelist:   /g/data/dg9/insar/Data61_VL/toolbox/PyRate/tests/test_data/cropA/ifg_30

# The DEM file used in the InSAR processing
demfile:       /g/data/dg9/insar/Data61_VL/toolbox/PyRate/tests/test_data/cropA/geotiffs/cropA_T005A_dem.tif

# The DEM header file from GAMMA (*.par) or ROI_PAC (*.rsc).
demHeaderFile: /g/data/dg9/insar/Data61_VL/toolbox/PyRate/tests/test_data/cropA/headers/cropA_20180106_VV_8rlks_eqa_dem.par

# File listing the pool of available header files (GAMMA: *mli.par, ROI_PAC: *.rsc)
hdrfilelist:   /g/data/dg9/insar/Data61_VL/toolbox/PyRate/tests/test_data/cropA/headers_13

# File listing the pool of available coherence files.
cohfilelist:   /g/data/dg9/insar/Data61_VL/toolbox/PyRate/tests/test_data/cropA/coherence_30

# File listing the pool of available baseline files (GAMMA).
basefilelist:  /g/data/dg9/insar/Data61_VL/toolbox/PyRate/tests/test_data/cropA/baseline_30

# Look-up table containing radar-coded row and column for lat/lon pixels (GAMMA)
ltfile:        /g/data/dg9/insar/Data61_VL/toolbox/PyRate/tests/test_data/cropA/geometry/20180106_VV_8rlks_eqa_to_rdc.lt

# Directory to write the outputs to
outdir:        {0}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PREPIFG parameters
#------------------------------------
# Input data format: ROI_PAC = 0, GAMMA = 1
processor:     1

# Coherence threshold value for masking, between 0 and 1
cohthresh:     0.05

# Multi-look/subsampling factor in east (x) and north (y) dimension
ifglksx:       1
ifglksy:       1

# Cropping options
# ifgcropopt: 1 = minimum extent 2 = maximum extent 3 = crop 4 = no cropping
# ifgxfirst,ifgyfirst: longitude (x) and latitude (y) of north-west corner
# ifgxlast,ifgylast: longitude (x) and latitude (y) of south-east corner
ifgcropopt:    1
ifgxfirst:     150.92
ifgyfirst:     -34.18
ifgxlast:      150.94
ifgylast:      -34.22

# No-data averaging threshold (0 = 0%; 1 = 100%)
noDataAveragingThreshold: 0.5

# The No-data value used in the interferogram files
noDataValue:   0.0

# Nan conversion flag. Set to 1 if missing No-data values are to be converted to NaN
nan_conversion: 1

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CORRECT parameters
#------------------------------------
# Reference pixel search options

# refx/y: Lon/Lat coordinate of reference pixel. If left blank then search for best pixel will be performed
# refnx/y: number of search grid points in x/y image dimensions
# refchipsize: size of the data window at each search grid point
# refminfrac: minimum fraction of valid (non-NaN) pixels in the data window
refx:          -1
refy:          -1
refnx:         5
refny:         5
refchipsize:   5
refminfrac:    0.8

#------------------------------------
# Reference phase correction method

# refest: 1 = median of the whole interferogram
# refest: 2 = median within the window surrounding the chosen reference pixel
refest:        1

#------------------------------------
# Orbital error correction

# orbfitmethod = 1: interferograms corrected independently; 2: network method
# orbfitdegrees: Degree of polynomial surface to fit (1 = planar; 2 = quadratic; 3 = part-cubic)
# orbfitlksx/y: additional multi-look factor for network orbital correction
orbfitmethod:  1
orbfitdegrees: 1
orbfitlksx:    1
orbfitlksy:    1

# phase closure params
# large_deviation_threshold_for_pixel # multiples of pi. For ex: 0.5 means pi/2, 1 means pi etc.
# avg_ifg_err_thr: ifgs with more than this fraction of pixels with error will be dropped
# loops_thr_ifg: pixel with phase unwrap error in at least this many loops
# phs_unw_err_thr: pixel with phase unwrap error in more than this many ifgs will be flagged
# max_loop_length: loops upto this many edges are considered for closure checks
# subtract_median: whether to subtract median during closure checks
# max_loops_in_ifg: loops with more than these many ifgs are discarded.
# max_loops_in_ifg: Ifg count must be met for all ifgs in the loop for loop to be discarded
large_dev_thr: 0.5
avg_ifg_err_thr: 0.07
loops_thr_ifg: 2
phs_unw_err_thr: 5
max_loop_length: 4
subtract_median: 1
max_loops_in_ifg: 2


#------------------------------------
# APS filter parameters

# tlpfcutoff: cutoff t0 for temporal high-pass Gaussian filter in days (int);
# tlpfpthr: valid pixel threshold;
# slpfcutoff: spatial low-pass Gaussian filter cutoff in km (greater than zero).
#             slpfcutoff=0 triggers cutoff estimation from exponential covariance function
tlpfcutoff:    30
tlpfpthr:      1
slpfcutoff:    1

#------------------------------------
# DEM error (residual topography) correction parameters

# de_pthr: valid observations threshold;
de_pthr:       20

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# TIMESERIES parameters
#------------------------------------

# tsmethod: Method for time series inversion (1 = Laplacian Smoothing; 2 = SVD)
# smorder: order of Laplacian smoothing operator (1 = first-order difference; 2 = second-order difference)
# smfactor: smoothing factor for Laplacian smoothing (value provided is converted as 10**smfactor)
# ts_pthr: valid observations threshold for time series inversion
tsmethod:      2
smorder:       2
smfactor:      -0.25
ts_pthr:       10

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# STACK parameters
#------------------------------------

# pthr: threshold for minimum number of ifg observations for each pixel
# nsig: threshold for iterative removal of observations
# maxsig: maximum sigma (std dev) used as an output masking threshold applied in Merge step. 0 = OFF.
pthr:          5
nsig:          3
maxsig:        100

# LOS Projection of timeseries and stack products
# converts slanted (native) LOS signals to either "vertical" or "horizontal", by multiplying by the sine or cosine of
# the incidence angle for each pixel.
# Three options, 0/1/2 - 0 = LOS (no conversion); 1 = pseudo-vertical; 2 = pseudo-horizontal.

los_projection: 0

# optionally supply rows and cols for tiles used during correct and merge step
rows:          4
cols:          3

largetifs:     0

[correct]
steps =
    orbfit
    refphase
    demerror
    phase_closure
    mst
    apscorrect
    maxvar
"""

def run_pyrate_cmd(pyrate_cmd, config_file, *args):
    """Run pyrate_cmd using config_file and any extra args."""
    cmd = ["pyrate", pyrate_cmd, "-f", config_file, *args]
    print("Running PyRate:", cmd)
    return subprocess.run(cmd, check=True, text=True)


# Create a temporary directory for the config and output files.
with TemporaryDirectory() as temp_output_dir:
    print("Created PyRate output directory", temp_output_dir)

    # Create the config file with interpolated values
    conf_file = os.path.join(temp_output_dir, "pyrate_job.conf")
    with open(conf_file, 'w') as f:
        f.write(CONF_FILE.format(temp_output_dir))

    # Run complete PyRate workflow from inside the pyrate toolbox directory
    # so the relative paths in the config files work.
    os.chdir("/g/data/dg9/insar/Data61_VL/toolbox/PyRate")
    try:
        run_pyrate_cmd("workflow", conf_file)
    except subprocess.CalledProcessError as ex:
        with open(os.path.join(temp_output_dir, "EXCEPTION.txt"), 'w') as f:
            f.write(f"PyRate job failed with exception.\nreturncode: {ex.returncode}\nfailed command: {ex.cmd}\n\nexception: {ex!r}")

    # Work around "cloud" bug in vgl by running in the output dir so we can use
    # relative filenames.
    os.chdir(temp_output_dir)
    # Zip directories for cloud upload
    for f in os.listdir(temp_output_dir):
        abs_f = os.path.join(temp_output_dir, f)
        if os.path.isdir(abs_f):
            fzip = f"{f}.zip"
            print("Zipping dir", f, "into", fzip, ".")
            with ZipFile(f"{f}.zip", 'a') as z:
                z.write(f)
    # Upload results, files only
    for f in os.listdir(temp_output_dir):
        abs_f = os.path.join(temp_output_dir, f)
        if not os.path.isdir(abs_f):
            # Work around vgl "cloud" bug by using relative filenames
            subprocess.run(["cloud", "upload", f, f])
            # subprocess.run(["cloud", "upload", f, abs_f])
