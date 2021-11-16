import os
import os.path
import subprocess
import sys
from tempfile import TemporaryDirectory

CONF_FILE = r"""
# PyRate configuration file for GAMMA-format interferograms
#
#------------------------------------
# input/output parameters

# Directory for the (unwrapped) interferograms.
obsdir:       ${obs_dir}

# File containing the list of interferograms to use.
ifgfilelist:  ${ifg_list}

# The DEM file used in the InSAR processing
demfile:      ${dem_file}

# The DEM header file from GAMMA (*.par) or ROI_PAC (*.rsc).
demHeaderFile: ${dem_header_file}

# GAMMA only: The directory containing GAMMA slc.par header files for all epochs
slcFileDir:   ${slc_file_dir}

# GAMMA only: File listing the pool of available slc.par header files
slcfilelist: ${slc_file_list}

# Directory containing the coherence files. If not provided, obsdir will be used.
cohfiledir: ${coh_file_dir}

# File listing the pool of available coherence files.
cohfilelist: ${coh_file_list}

# Directory to write the outputs to
outdir:       {0}

# InSAR processing software: ROI_PAC = 0, GAMMA = 1
processor:    ${processor}

# No data averaging threshold for prepifg
noDataAveragingThreshold: ${nodata_avg_thr}

# The no data value in the interferograms
noDataValue:  ${nodata_value}

# Nan conversion flag. Set to 1 if missing (0) phase values are converted to nan
nan_conversion: ${nan_conv}

#-----------------------------------
# Multi-threading parameters: used by linrate/timeseries/prepifg
# gamma prepifg runs in parallel on a single machine if parallel != 0
# parallel = 1, linrate/timeseries computation is done in parallel by rows
# parallel = 2, linrate/timeseries computation is done in parallel for each pixel
# parallel = 0, linrate/timeseries computation is done in serial pixel by pixel
parallel:  ${parallel}
processes: ${processes}

#------------------------------------
# Coherence masking options: used by process
# cohmask: 1 = ON, 0 = OFF
# cohthresh: coherence threshold value, between 0 and 1
cohmask:   ${cohmask}
cohthresh:  ${cohthresh}

#------------------------------------
# Interferogram multi-look and crop options
# ifgcropopt: 1 = minimum 2 = maximum 3 = customise 4 = all ifms already same size
# ifglksx/y: multi-look/subsampling factor in east and north direction respectively
# ifgxfirst,ifgyfirst: x,y of top-left corner
# ifgxlast,ifgylast: x,y of bottom-right corner
ifgcropopt:   ${ifgcropopt}
ifglksx:      ${ifglksx}
ifglksy:      ${ifglksy}
ifgxfirst:    ${ifgxfirst}
ifgxlast:     ${ifgxlast}
ifgyfirst:    ${ifgyfirst}
ifgylast:     ${ifgylast}

#------------------------------------
# Reference pixel search options
# refx/y: Lon/Lat coordinate of reference pixel. If left blank then search for best pixel will be performed
# refnx/y: number of search grid points in x/y direction
# refchipsize: chip size of the data window at each search grid point
# refminfrac: minimum fraction of valid (non-NaN) pixels in the data window
refx:          ${refx}
refy:          ${refy}
refnx:         ${refnx}
refny:         ${refny}
refchipsize:   ${refchipsize}
refminfrac:    ${refminfrac}

#------------------------------------
# Reference phase calculation method
# refest: 1 = median of the whole interferogram
# refest: 2 = median within the window surrounding the chosen reference pixel
refest:        ${refest}

#------------------------------------
# Orbital error correction
# orbfit: ON = 1, OFF = 0
# orbfitmethod = 1: independent method; 2: network method
# orbfitdegrees: Degree of polynomial surface to fit (1 = planar; 2 = quadratic; 3 = part-cubic)
# orbfitlksx/y: additional multi-look factor for orbital correction
orbfit:        ${orbfit}
orbfitmethod:  ${orbfitmethod}
orbfitdegrees: ${orbfitdegrees}
orbfitlksx:    ${orbfitlksx}
orbfitlksy:    ${orbfitlksy}

#------------------------------------
# APS correction using spatio-temporal filter
# apsest: ON = 1, OFF = 0
# Spatial low-pass filter parameters
# slpfmethod: filter method (1: butterworth; 2: gaussian)
# slpfcutoff: cutoff d0 (greater than zero) in km for both butterworth and gaussian filters
# slpforder: order n for butterworth filter (default 1)
# slpnanfill: 1 for interpolation, 0 for zero fill
# slpnanfill_method: linear, nearest, cubic; only used when slpnanfill=1
# Temporal low-pass filter parameters
# tlpfmethod: 1 = Gaussian, 2 = Triangular, 3 = Mean filter
# tlpfcutoff: cutoff t0 for gaussian filter in year;
# tlpfpthr: valid pixel threshold;
apsest:         ${apsest}
slpfmethod:     ${slpfmethod}
slpfcutoff:     ${slpfcutoff}
slpforder:      ${slpforder}
slpnanfill:     ${slpnanfill}
slpnanfill_method:  ${slpnanfill_method}
tlpfmethod:   ${tlpfmethod}
tlpfcutoff:   ${tlpfcutoff}
tlpfpthr:     ${tlpfpthr}

#------------------------------------
# Time Series Calculation
# tscal: YES = 1, NO = 0
# tsmethod: Method for time series inversion (1 = Laplacian Smoothing; 2 = SVD)
# smorder: order of Laplacian smoothing operator (1 =  first-order difference; 2 = second-order difference)
# smfactor: smoothing factor for Laplacian smoothing (value provided is converted as 10**smfactor)
# ts_pthr: valid observations threshold for time series inversion
tscal:         ${tscal}
tsmethod:      ${tsmethod}
smorder:       ${smorder}
smfactor:      ${smfactor}
ts_pthr:       ${ts_pthr}

#------------------------------------
# Stacking calculation
# pthr: minimum number of coherent ifg connections for each pixel
# nsig: n-sigma used as residuals threshold for iterativel least squares stacking
# maxsig: maximum residual used as a threshold for values in the rate map
nsig:          ${nsig}
pthr:          ${pthr}
maxsig:        ${maxsig}
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

    # Run PyRate
    try:
        run_pyrate_cmd("conv2tif", conf_file)
        run_pyrate_cmd("prepifg", conf_file)
        run_pyrate_cmd("process", conf_file)
        run_pyrate_cmd("merge", conf_file)
    except subprocess.CalledProcessError as ex:
        with open(os.path.join(temp_output_dir, "EXCEPTION.txt"), 'w') as f:
            f.write(f"PyRate job failed with exception.\nreturncode: {ex.returncode}\nfailed command: {ex.cmd}\n\nexception: {ex!r}")

    # Work around "cloud" bug in vgl by running in the output dir so we can use
    # relative filenames.
    os.chdir(temp_output_dir)
    # Upload results, files only
    for f in os.listdir(temp_output_dir):
        abs_f = os.path.join(temp_output_dir, f)
        if not os.path.isdir(abs_f):
            # Work around vgl "cloud" bug by using relative filenames
            subprocess.run(["cloud", "upload", f, f])
            # subprocess.run(["cloud", "upload", f, abs_f])
