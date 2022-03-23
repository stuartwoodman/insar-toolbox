import os
import os.path
import subprocess
import sys
from glob import glob
from tempfile import TemporaryDirectory
import xml.etree.ElementTree as ET
from datetime import date
import re
from typing import List
from osgeo import ogr

# ================
# VIC and NSW 2021
# ================
INSAR_ARD_DIR = "/g/data/dz56/insar_final_processing/"
INSAR_BASE_VARIANT = "_VV_8rlks"
INSAR_ARD_VARIANT = "_VV_8rlks_geo"
INSAR_COH_VARIANT = "_VV_8rlks_flat_geo"

# ========
# VICTORIA
# ========
# INSAR_ARD_DIR = "/g/data/dz56/INSAR_ARD/VV/INSAR_ANALYSIS/VICTORIA/S1/GAMMA/"
# INSAR_BASE_VARIANT = "_VV_4rlks"
# INSAR_ARD_VARIANT = "_VV_4rlks_eqa"
# INSAR_COH_VARIANT = "_VV_4rlks_flat_eqa"

INSAR_INTERVAL_DIR_RE = re.compile("(\d\d\d\d)(\d\d)(\d\d)-(\d\d\d\d)(\d\d)(\d\d)")
INSAR_DEM_RE = re.compile("(\d\d\d\d)(\d\d)(\d\d).*")
INSAR_FRAME_RE = re.compile("T(\d+)D_F(\d+)S_(.+)")

ISO_DATE_RE = re.compile("(\d\d\d\d)-(\d\d?)-(\d\d?)")

def parse_iso_date(iso_date):
    """Return the date object parsed from iso_date or None.

    Return None if iso_date is None or the empty string. A non-empty string that
    fails to parse will raise an error.

    """
    if iso_date:
        match = ISO_DATE_RE.match(iso_date)
        if match:
            return date(*map(int, match.group(1,2,3)))
        print("Configured date not a valid ISO8601 string:", iso_date)
    return None

# TODO: Use start/end dates to select interferograms, but how to resolve date
# ranges in input files where there are overlaps? E.g for the following three
# ranges, assuming start <= 20180606 and end >= 20180630, should we use [1,3],
# [2], [1,2,3]?
#
# 1. 20180606,20180618
# 2. 20180606,20180630
# 3. 20180618,20180630
#
# Answer from MattG: All of them. Use all ifgs with intervals that are completely within the specified interval.
#
# Dates will be parsed as ISO8601 strings, and assumed to match the timezone of
# the date strings used in the interferogram path/filename structures.
#
START_DATE_ISO = "${startdate}"
END_DATE_ISO = "${enddate}"
start_date = parse_iso_date(START_DATE_ISO)
end_date = parse_iso_date(END_DATE_ISO)
print("Configured start date:", start_date.isoformat() if start_date else "None, will use all available data.")
print("Configured end date:", end_date.isoformat() if end_date else "None, will use all available data.")

# Tile dataset to grab pyrate tile data from.
insar_tiles = "${insar_tiles}"

# Bounding box to crop processing in geojson bbox format "west lon, south lat,
# east lon, north lat". Values can be separated by commas and/or whitespace.
BBOX_RE = re.compile("\[?(-?\d+(?:\.\d+)?)[\s,]+(-?\d+(?:\.\d+)?)[\s,]+(-?\d+(?:\.\d+)?)[\s,]+(-?\d+(?:\.\d+)?)\]?")
crop_bbox = "${crop_bbox}"

CONF_FILE = r"""
# PyRate configuration file for GAMMA-format interferograms
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Optional ON/OFF switches - ON = 1; OFF = 0

# Coherence masking (PREPIFG)
cohmask:   ${cohmask}

# Orbital error correction (CORRECT)
orbfit:    ${orbfit}

# APS correction using spatio-temporal filter (CORRECT)
apsest:    ${apsest}

# DEM error (residual topography) correction (CORRECT)
demerror:      ${demerror}

# Phase Closure correction (CORRECT)
phase_closure: 1

# Time series calculation (PROCESS)
tscal:     ${tscal}

# Optional save of numpy array files for output products (MERGE)
savenpy:       ${savenpy}

# Optional save of incremental time series products (TIMESERIES/MERGE)
savetsincr:    ${savetsincr}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Integer parameters

# LOS Projection of output products converts slanted (native) LOS signals
# to either "pseudo-vertical" or "pseudo-horizontal", by dividing by the
# cosine or sine of the incidence angle for each pixel, respectively.
# los_projection: 0 = LOS (no conversion); 1 = pseudo-vertical; 2 = pseudo-horizontal.
los_projection: ${los_projection}

# Number of sigma to report velocity error. Positive integer. Default: 2 (TIMESERIES/STACK)
velerror_nsig: 2

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Multi-threading parameters used by stacking/timeseries/prepifg
# gamma prepifg runs in parallel on a single machine if parallel = 1
# parallel: 1 = parallel, 0 = serial
parallel:  0
# number of processes
processes: 1

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Input/Output file locations
#
# File containing the list of interferograms to use.
ifgfilelist:   {ifgfilelist}

# The DEM file used in the InSAR processing
demfile:       {demfile}

# The DEM header file from GAMMA (*.par) or ROI_PAC (*.rsc).
demHeaderFile: {demheaderfile}

# File listing the pool of available header files (GAMMA: *slc.par, ROI_PAC: *.rsc)
hdrfilelist:   {hdrfilelist}

# File listing the pool of available coherence files.
cohfilelist:   {cohfilelist}

# File listing the pool of available baseline files (GAMMA).
#basefilelist:  tests/test_data/cropA/baseline_30
basefilelist:  {basefilelist}

# Look-up table containing radar-coded row and column for lat/lon pixels (GAMMA)
#ltfile:        tests/test_data/cropA/geometry/20180106_VV_8rlks_eqa_to_rdc.lt
ltfile:        {ltfile}

# Directory to write the outputs to
outdir:        {outdir}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PREPIFG parameters
#------------------------------------
# Input data format: ROI_PAC = 0, GAMMA = 1
processor:    1

# Coherence threshold value for masking, between 0 and 1
cohthresh:  ${cohthresh}

# Multi-look/subsampling factor in east (x) and north (y) dimension
ifglksx:      ${ifglksx}
ifglksy:      ${ifglksy}

# Cropping options
# ifgcropopt: 1 = minimum extent 2 = maximum extent 3 = crop 4 = no cropping
# ifgxfirst,ifgyfirst: longitude (x) and latitude (y) of north-west corner
# ifgxlast,ifgylast: longitude (x) and latitude (y) of south-east corner
ifgcropopt:   3
{ifgxfirst}
{ifgyfirst}
{ifgxlast}
{ifgylast}

# No-data averaging threshold (0 = 0%; 1 = 100%)
noDataAveragingThreshold: 0.5

# The No-data value used in the interferogram files
noDataValue:  0.0

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
refx:          ${refx}
refy:          ${refy}
refnx:         ${refnx}
refny:         ${refny}
refchipsize:   ${refchipsize}
refminfrac:    ${refminfrac}

#------------------------------------
# Reference phase correction method

# refest: 1 = median of the whole interferogram
# refest: 2 = median within the window surrounding the chosen reference pixel
refest:        ${refest}

#------------------------------------
# Orbital error correction

# orbfitmethod = 1: interferograms corrected independently; 2: network method
# orbfitdegrees: Degree of polynomial surface to fit (1 = planar; 2 = quadratic; 3 = part-cubic)
# orbfitlksx/y: additional multi-look factor for orbital correction
orbfitmethod:  ${orbfitmethod}
orbfitdegrees: ${orbfitdegrees}
orbfitlksx:    10
orbfitlksy:    10

#------------------------------------
# Phase closure correction parameters

# large_deviation_threshold_for_pixel # multiples of pi. For ex: 0.5 means pi/2, 1 means pi etc.
# avg_ifg_err_thr: ifgs with more than this fraction of pixels with error will be dropped
# loops_thr_ifg: pixel with phase unwrap error in at least this many loops
# phs_unw_err_thr: pixel with phase unwrap error in more than this many ifgs will be flagged
# max_loop_length: loops upto this many edges are considered for closure checks
# subtract_median: whether to subtract median during closure checks
# max_loops_in_ifg: loops with more than these many ifgs are discarded.
# max_loops_in_ifg: Ifg count must be met for all ifgs in the loop for loop to be discarded
large_dev_thr:    0.5
avg_ifg_err_thr:  0.07
loops_thr_ifg:    2
phs_unw_err_thr:  5
max_loop_length:  4
subtract_median:  1
max_loops_in_ifg: 2


#------------------------------------
# APS spatial low-pass filter parameters

# slpfmethod: Spatial low-pass filter method (1: butterworth; 2: gaussian)
# slpfcutoff: cutoff d0 (greater than zero) in km for both butterworth and gaussian filters
# slpforder: order n for butterworth filter (default 1)
# slpnanfill: 1 for interpolation, 0 for zero fill
# slpnanfill_method: linear, nearest, cubic; only used when slpnanfill=1
slpfmethod:     2
slpfcutoff:     ${slpfcutoff}
slpforder:      1
slpnanfill:     1
slpnanfill_method:  cubic

#------------------------------------
# APS temporal low-pass filter parameters

# tlpfmethod: 1 = Gaussian, 2 = Triangular, 3 = Mean filter
# tlpfcutoff: cutoff t0 for gaussian filter in year;
# tlpfpthr: valid pixel threshold;
tlpfmethod:   1
tlpfcutoff:   ${tlpfcutoff}
tlpfpthr:     1

#------------------------------------
# Time Series Calculation parameters

# tsmethod: Method for time series inversion (1 = Laplacian Smoothing; 2 = SVD)
# smorder: order of Laplacian smoothing operator (1 = first-order difference; 2 = second-order difference)
# smfactor: smoothing factor for Laplacian smoothing (value provided is converted as 10**smfactor)
# ts_pthr: valid observations threshold for time series inversion
tsmethod:      2
smorder:       2
smfactor:     -0.25
ts_pthr:       ${ts_pthr}

#------------------------------------
# Stacking calculation parameters

# pthr: minimum number of coherent ifg connections for each pixel
# nsig: n-sigma used as residuals threshold for iterative least squares stacking
# maxsig: maximum residual used as a threshold for values in the rate map
pthr:          ${pthr}
nsig:          ${nsig}
maxsig:        ${maxsig}

#------------------------------------
# DEM error (residual topography) correction parameters

# de_pthr: valid observations threshold;
de_pthr:       20
"""

# Namespaces for wfs/gml/etc
NS_WFS1 = {
    "wfs": "http://www.opengis.net/wfs",
    "gml": "http://www.opengis.net/gml",
    "insar": "http://csiro.au/insar"
}
PATH_WFS1 = "./gml:featureMembers//insar:insar_tiles"
# TODO: look for insar tile schema rather than element name

NS_WFS2 = {
    "wfs": "http://www.opengis.net/wfs/2.0",
    "gml": "http://www.opengis.net/gml/3.2",
    "insar": "http://csiro.au/insar"
}
PATH_WFS2 = "./wfs:member//insar:insar_tiles"
# TODO: look for insar tile schema rather than element name

PATH_GEOMETRY = "./insar:wkb_geometry/gml:Polygon/gml:exterior/gml:LinearRing/gml:posList"


class PyrateException(Exception):
    pass


class InsarInterval(object):
    """Wrapper for an interferogram interval. """
    def __init__(self, name=None, path=None, start=None, end=None,
                 unw=None, tif=None):
        self.name = name
        self.path = path
        self.start = start
        self.end = end
        self.unw = unw
        self.tif = tif


class InsarTile(object):
    """Wrapper for an INSAR tile feature."""

    def __init__(self, tile, ns):
        self.tile = tile
        self.ns = ns
        self._intervals = []
        self.gmlid = tile.get(f"{{{ns['gml']}}}id")

    @property
    def relorb(self):
        """Return the relative orbit."""
        return self._prop("insar:track")

    @property
    def track(self):
        """Return the relative orbit."""
        return self.relorb

    @property
    def frame(self):
        """Return the frame. """
        return self._prop("insar:frame_id")

    @property
    def tiledir(self):
        """Return the name of the tile directory."""
        return self._prop("insar:tile_name")

    def scan_ard(self, ard_dir, variant, coh_variant, base_variant, start=None, end=None):
        """Scan INSAR files for this tile and return number of interferograms.

        If a start and/or end date are specified then only include the
        interferograms that fall within the specified date range.

        """
        # Check we have an interferograms directory
        ifgdir = os.path.join(ard_dir, "INT")
        if not os.path.isdir(ifgdir):
            raise PyrateException(f"Interferograms directory not found: {ifgdir}")

        # Find all interferogram directories in the ard_dir, and sort by date
        self._intervals = []
        with os.scandir(ifgdir) as it:
            for entry in it:
                if entry.is_dir():
                    self.load_interval(entry, variant, coh_variant, base_variant, start, end)
        self._intervals = sorted(self._intervals, key=lambda x: x.start)

        # Find the DEM files for this tile
        demdir = os.path.join(ard_dir, "DEM")
        if not os.path.isdir(demdir):
            print("DEM directory not found in standard location.")
        else:
            print("Standard DEM directory found for tile:", demdir)
            # Find latest dem file, looking for a .dem.tif version first then
            # raw .dem files.
            latest = None
            dem_f = None
            dem_h = None
            rdc_lt = None
            globs = [*glob(os.path.join(demdir, f"*{variant}.dem.tif")),
                     *glob(os.path.join(demdir, f"*{variant}.dem"))]
            for f in globs:
                demdate = INSAR_DEM_RE.match(os.path.basename(f))
                if latest is None or demdate > latest:
                    latest = demdate
                    dem_f = f
            if not dem_f:
                print("No DEM file found in DEM directory.")
            else:
                print("Found standard DEM file:", dem_f)
                self.dem_file = dem_f
                # DEM header can be *.dem.par or *_dem.par, where the * is the
                # same as the stem from the DEM file, but the '.' or '_' before
                # 'dem.par' may be different to the pattern used in the DEM
                # file!
                #
                # Find the stem from dem_f, strip the trailing '.' or '_', and
                # append '[._]dem.par' to find the header file.
                stem = dem_f.rpartition(".dem")[0]
                for c in "._":
                    tryh = f"{stem}{c}dem.par"
                    if os.path.isfile(tryh):
                        dem_h = tryh
                        break
                if dem_h:
                    print("Found matching standard DEM header file:", dem_h)
                    self.dem_header = dem_h
                else:
                    print("No DEM header file found to match DEM file")
                # Repeat the search for the "<stem>_to_rdc.lt" lookup table file.
                trylt = f"{stem}_to_rdc.lt"
                if os.path.isfile(trylt):
                    print("Found matching standard lookup table file:", trylt)
                    self.lt_file = trylt
                else:
                    print("No matching lookup table file found to match DEM file.")

        # Find the SLC header files
        slcdir = os.path.join(ard_dir, "SLC")
        if not os.path.isdir(slcdir):
            print("SLC headers directory not found in standard location.")
        else:
            print("Found standard SLC directory for tile:", slcdir)
            self.slcs = glob(os.path.join(slcdir, "*", "*_VV.slc.par"))

    def load_interval(self, entry, variant, coh_variant, base_variant, start_date=None, end_date=None):
        """Load the interval contained in DirEntry entry.

        If we have a configured start and/or end date then the interval will
        only be loaded if it falls entirely within the specified range.

        """
        # Parse directory name for start/end of the interval in the directory.
        match = INSAR_INTERVAL_DIR_RE.match(entry.name)
        if match:
            start = date(*map(int, match.group(1,2,3)))
            end = date(*map(int, match.group(4,5,6)))

            # If we have a temporal constraint and this interval falls outside
            # it then do not load it.
            if ((start_date and start < start_date) or
                (end_date and end > end_date)):
                print("Ignoring interval directory outside configured date range:", entry.name)
                return

            # This interval is OK so continue to load it.
            print("Loading interval directory:", entry.name)
            interval = InsarInterval(entry.name, entry.path, start, end)
            unw = os.path.join(entry.path, f"{entry.name}{variant}.unw")
            if os.path.isfile(unw):
                interval.unw = unw
            tif = os.path.join(entry.path, f"{entry.name}{variant}.unw.tif")
            if os.path.isfile(tif):
                interval.tif = tif
            else:
                # Try the other tif filename style
                tif = os.path.join(entry.path, f"{entry.name}{variant}_unw.tif")
                if os.path.isfile(tif):
                    interval.tif = tif
            if interval.unw or interval.tif:
                # Found an interferogram file.
                # Find the matching coherence file
                coh = os.path.join(entry.path, f"{entry.name}{coh_variant}.coh.tif")
                if os.path.isfile(coh):
                    interval.coh = coh
                else:
                    # Try underscore style
                    coh = os.path.join(entry.path, f"{entry.name}{coh_variant}_coh.tif")
                    if os.path.isfile(coh):
                        interval.coh = coh
                    else:
                        print("Missing coherence file", coh, "for interferogram:", interval.tif if interval.tif else interval.unw)
                # Find the matching baseline file
                basefile = os.path.join(entry.path, f"{entry.name}{base_variant}_base.par")
                if os.path.isfile(basefile):
                    interval.basefile = basefile
                else:
                    print("Missing baseline file", basefile, "for interferogram:", interval.tif if interval.tif else interval.unw)
                self._intervals.append(interval)
            else:
                print(f"No unwrapped interferogram or tiff found in {entry.path}")
        else:
            print("Directory", entry.name, "Does not match INSAR interval pattern.")

    def intervals(self, start=None, end=None):
        """Return list of intervals for this tile between start and end.

        If start is None return from the earliest available.
        If end is None return to the latest available.

        An interval is returned if it is falls within [start, end], not if it
        overlaps before or after the specified temporal range.

        """
        if self._intervals:
            if start is None:
                start = self._intervals[0].start
            if end is None:
                end = self._intervals[-1].end
        return [i for i in self._intervals if i.start >= start and i.end <= end]

    def _prop(self, element):
        """Return the content of element in tile."""
        elem = self.tile.find(element, namespaces=self.ns)
        if elem is not None:
            return elem.text
        return None


class InsarTileFeatures(object):
    """Parse a feature collection of INSAR tiles."""

    def load(self, tiles_path):
        """Load and parse the collection from file at tiles_path."""
        self.tree = ET.parse(tiles_path)
        self.root = self.tree.getroot()
        if self.root.tag == "{http://www.opengis.net/wfs}FeatureCollection":
            print("Found WFS v1 FeatureCollection at root.")
            self.ns = NS_WFS1
            self.PATH = PATH_WFS1
        elif self.root.tag == "{http://www.opengis.net/wfs/2.0}FeatureCollection":
            print("Found WFS v2 FeatureCollection at root.")
            self.ns = NS_WFS2
            self.PATH = PATH_WFS2
        else:
            raise PyrateException("Not a WFS FeatureCollection, giving up.")
        print("Namespace map:", self.ns)
        print("Tiles path:", self.PATH)

    def tiles(self):
        """Return list of all tiles in the collection."""
        if not self.root:
            raise Exception("PyrateTiles.tiles() called before load().")

        tiles = self.root.findall(self.PATH, self.ns)
        if tiles:
            return [InsarTile(tile, self.ns) for tile in tiles]
        raise PyrateException("No tile found in input dataset.")

    def first(self):
        """Return the first tile in the collection."""
        if not self.root:
            raise Exception("PyrateTiles.first() called before load().")

        tile = self.root.find(self.PATH, self.ns)
        if tile:
            return InsarTile(tile, self.ns)
        raise PyrateException("No tile found in input dataset.")

    def bbox_to_wkt_text(self, bbox):
        """Return a WKT Polygon string from the crop bounding box"""
        point_list = bbox.split(', ')
        wkt = "Polygon(({x1} {y1},{x2} {y1},{x2} {y2},{x1} {y2},{x1} {y1}))".format(
            x1=point_list[0], y1=point_list[1], x2=point_list[2], y2=point_list[3]
        )
        return wkt

    def geom_to_wkt_text(self, geom):
        """Return a WKT Polygon string from the tile geometry"""
        point_list = geom.split(' ')
        wkt = "Polygon(("
        for i in range(0, len(point_list), 2):
            wkt += point_list[i] + " " + point_list[i+1]
            if i < len(point_list) - 2:
                wkt += ","
        wkt += "))"
        return wkt

    def best(self):
        """Return the tile in the collection completely constrained by crop_bbox."""
        if not self.root:
            raise Exception("PyrateTiles.best() called before load().")
        crop_wkt = self.bbox_to_wkt_text(crop_bbox)
        crop_poly = ogr.CreateGeometryFromWkt(crop_wkt)
        tiles = self.root.findall(self.PATH, self.ns)
        if tiles:
            for t in tiles:
                geom = t.find(PATH_GEOMETRY, self.ns)
                tile_wkt = self.geom_to_wkt_text(geom.text)
                tile_poly = ogr.CreateGeometryFromWkt(tile_wkt)
                if tile_poly.Contains(crop_poly):
                    return InsarTile(t, self.ns)
        raise PyrateException("No tile could be found that fully encloses the crop bounding box.")


# Return epoch strings in pyrate filename/path.
#
# Copied from https://github.com/GeoscienceAustralia/PyRate/blob/c3054c569b4b1827b61326741277606bc47e8463/pyrate/core/shared.py#L1373
def extract_epochs_from_filename(filename_with_epochs: str) -> List[str]:
    src_epochs = re.findall(r"(\d{8})", str(filename_with_epochs))
    if not len(src_epochs) > 0:
        src_epochs = re.findall(r"(\d{6})", str(filename_with_epochs))
    return src_epochs


def find_epoch(f):
    """Return the first epoch string extracted from filename f."""
    if f:
        epochs = extract_epochs_from_filename(f)
        if len(epochs) > 0:
            return epochs[0]
    return None


dataset = InsarTileFeatures()
try:
    # Load tile features from input dataset
    dataset.load(insar_tiles)

    # Grab the first tile
    # TODO: Support processing multiple tiles
    tile = dataset.best()
    tile_id = tile.gmlid
    print("Found tile", tile_id)
except Exception as ex:
    print("Loading INSAR tile features failed:", ex)
    sys.exit(1)

tile_dir = os.path.join(INSAR_ARD_DIR, tile.tiledir)
if not os.path.isdir(tile_dir):
    print("Interferograms directory from tile info not found")
    print(tile_dir, " is not a directory.")
    sys.exit(1)
print("INSAR ARD data from", tile_dir)


def run_bad_interferogram_detector(input_file, output_file):
    """Run bad interferogram detector """
    cmd = ["python3", "/g/data/dg9/rlt118/bad_int_detector/stack_filter_BID.py", input_file, output_file]
    print("Running bad interferogram detector")
    return subprocess.run(cmd, check=True, text=True)


def run_pyrate_cmd(pyrate_cmd, config_file, *args):
    """Run pyrate_cmd using config_file and any extra args."""
    cmd = ["pyrate", pyrate_cmd, "-f", config_file, *args]
    print("Running PyRate:", cmd)
    return subprocess.run(cmd, check=True, text=True)


# Create a temporary directory for the config and output files.
with TemporaryDirectory() as temp_output_dir:
    print("Created PyRate output directory", temp_output_dir)

    # Scan the data directory for the tile
    tile.scan_ard(tile_dir,
                  INSAR_ARD_VARIANT,
                  INSAR_COH_VARIANT,
                  INSAR_BASE_VARIANT,
                  start_date,
                  end_date)
    intervals = tile.intervals()
    has_tifs = all(i.tif for i in intervals)
    if has_tifs:
        print("Found tifs for all interferograms, will skip conv2tif")

    # Generate the ifg list
    ifgfilelist = os.path.join(temp_output_dir, "interferograms.list")
    cohfile = os.path.join(temp_output_dir, "coh_list.txt")
    basefilelist = os.path.join(temp_output_dir, "baselinefiles.list")
    with open(ifgfilelist, 'w') as f:
        with open(cohfile, 'w') as g:
            with open(basefilelist, 'w') as h:
                for interval in intervals:
                    f.write(interval.tif if has_tifs else interval.unw)
                    f.write('\n')
                    g.write(interval.coh)
                    g.write('\n')
                    h.write(interval.basefile)
                    h.write('\n')

    # Run bad interferogram detector
    try:
        ifgcleanfilelist = os.path.join(temp_output_dir, "ifgs-bid.list")
        run_bad_interferogram_detector(ifgfilelist, ifgcleanfilelist)
    except subprocess.CalledProcessError as ex:
        print(f"PyRate bad interferogram detection failed with exception.\nreturncode: {ex.returncode}\nfailed command: {ex.cmd}\n\nexception: {ex!r}")

    # Find the DEM and its header file
    if not tile.dem_file or not tile.dem_header:
        print("No DEM file or header found in standard location.")
        print("I can't proceed without a DEM file and header specified manually.")
        sys.exit(1)

    # Write the header files list and coherence files list.
    if not tile.slcs:
        print("No SLC header files for tile.")
        print("I can't proceed without a list of header files.")
        sys.exit(1)
    hdrfile = os.path.join(temp_output_dir, "hdr_list.txt")
    with open(hdrfile, 'w') as f:
        # NB - See gamma_header() in pyrate/core/gamma.py for use of this header
        # file. It appears to assume the headers are sorted by epoch in
        # ascending date order.
        for slc in sorted(tile.slcs, key=find_epoch):
            f.write(slc)
            f.write('\n')

    config = dict(ifgfilelist=ifgcleanfilelist,
                  demfile=tile.dem_file,
                  demheaderfile=tile.dem_header,
                  ltfile=tile.lt_file,
                  basefilelist=basefilelist,
                  hdrfilelist=hdrfile,
                  cohfilelist=cohfile,
                  outdir=temp_output_dir,
                  ifgxfirst='',
                  ifgyfirst='',
                  ifgxlast='',
                  ifgylast='')

    # Parse the ifg crop lat/lon values from the crop_bbox string. NB pyrate
    # expects ifgyfirst to be the northerly latitude (i.e. they want TL-BR
    # ordering)
    match = BBOX_RE.match(crop_bbox.strip())
    if match:
        config['ifgxfirst'] = f"ifgxfirst:    {match.group(1)}"
        config['ifgyfirst'] = f"ifgyfirst:    {match.group(4)}"
        config['ifgxlast'] = f"ifgxlast:     {match.group(3)}"
        config['ifgylast'] = f"ifgylast:     {match.group(2)}"

    # Create the config file with interpolated values
    conf_file = os.path.join(temp_output_dir, "pyrate_job.conf")
    with open(conf_file, 'w') as f:
        f.write(CONF_FILE.format(**config))

    # Run PyRate if we have a hopefully valid config
    if all(v is not None for v in config.values()):
        try:
            # Run the full workflow
            # TODO allow user to configure the steps they want to run
            #run_pyrate_cmd("workflow", conf_file)
            if not has_tifs:
                run_pyrate_cmd("conv2tif", conf_file)
            run_pyrate_cmd("prepifg", conf_file)
            run_pyrate_cmd("correct", conf_file)
            run_pyrate_cmd("timeseries", conf_file)
            #run_pyrate_cmd("process", conf_file)
            run_pyrate_cmd("merge", conf_file)
        except subprocess.CalledProcessError as ex:
            with open(os.path.join(temp_output_dir, "EXCEPTION.txt"), 'w') as f:
                f.write(f"PyRate job failed with exception.\nreturncode: {ex.returncode}\nfailed command: {ex.cmd}\n\nexception: {ex!r}")

    # Upload results, files only
    # # Work around "cloud" bug in vgl by running in the output dir so we can use
    # # relative filenames.
    # os.chdir(temp_output_dir)

    # Upload utility
    # TODO: Pull out into VGL python library, replacing our old
    # cloud.sh script and the version from nci-util.sh
    def upload_results(d, prefix=None):
        """Upload results from directory d, with optional prefix.

        A result file is uploaded under its filename, possibly prefixed as
        specified.

        All the contents of a result directory will be uploaded with the
        directory name as a prefix. It will recurse into subdirectories as
        required.

        """
        if prefix is None:
            prefix = ""
        subprocess.run(["bash", "-c", f"source nci-util.sh; cloud upload output {d}"])
		
	
    json_path = os.path.join(temp_output_dir, "json")
    cmd = ["time_series_to_json.py", "-i", temp_output_dir, "-o", json_path, "-z"]
    subprocess.run(cmd, check=True, text=True)

    upload_results(temp_output_dir)
