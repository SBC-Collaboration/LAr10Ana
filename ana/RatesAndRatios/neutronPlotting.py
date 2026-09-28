## Variant of neutronExcluded.py: identical pipeline -- same livetime accounting, same combined_region()
## region assignment, same background dome-cut symmetry -- except build_groups()'s normalizationFactor
## goes back to being FIT to the measured rates, livetime-weighted per pipeline, instead of derived from
## the source activity alone. This is the normalization neutronExcluded.py used before it switched to
## activityNormalizationFactor; see that file's header for the activity-based version and why it was
## introduced (so a sim/data disagreement is a real disagreement rather than something a fit absorbed).
##
## The fit, done once per pipeline inside build_groups() after every group's raw sim counts are known:
##   ratio_i = sum(observed rate_i) / sum(raw sim count_i)      -- one ratio per group
##   normalizationFactor = sum(liveTime_i * ratio_i) / sum(liveTime_i)   -- livetime-weighted mean
## Groups with more exposure (a more reliable data/sim ratio) pull the shared factor harder than short,
## noisy ones. Run this alongside neutronExcluded.py and compare combinedMultiplicityFinal.pdf between
## them to see what the fit is absorbing that the activity-based number doesn't.
##
## Quick A/B comparison copy, not meant to replace neutronExcluded.py -- writes to its own PLOTS_ROOT
## ("plotsFittedNorm/") so it never overwrites the main pipeline's figures/reports. OUTPUT_ROOT
## ("output/") is left shared: write_final_rates()'s tables never reference normalizationFactor at all
## (they are raw event-count-over-exposure, uncut, with no simulation involved), so they come out
## byte-identical between the two scripts and there is no reason to duplicate them.
##
## dome exclusion: same rate/threshold-fit pipeline as neutronPlotting.py, but the sim and real sides both drop dome-region hits, gated by excludedRegions below
## On the real side the cut drops the BUBBLE and keeps the exposure's live time -- see bin_multiplicities().
import atexit
import glob
import os
import sys
from datetime import date

import matplotlib.pyplot as plt
import numpy as np
from sbcbinaryformat import Streamer
import SeitzModel as sm

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(SCRIPT_DIR, "neutronSim"))
CF_SIM_DIR = "/nashome/o/ochiarin/Documents/neutronSim"
sys.path.insert(0, CF_SIM_DIR)
from cfconfBThresholds import get_event_coordinates, MULTIPLICITY_CUT
from cfconfBThresholds import get_multiplicity_counts as cfconf_get_multiplicity_counts

## config variables
HANDSCAN_DIR = "/exp/e961/data/SBC-25-handscan/"
RECON_DIR = "/exp/e961/data/SBC-25-recon/v1.0.0/"

DOME_Z_THRESHOLD_CM = -3

# the Seitz threshold is evaluated at the measured PT1101 pressure rather than the nominal pset. Two steps:
# a linear calibration putting the raw transducer reading on the LAr-pressure scale at the transducer
# (p = raw * GAIN - OFFSET), then a fixed hydrostatic drop down to the target volume itself
PT1101_GAIN = 1.0388
PT1101_OFFSET = 0.589
PT1101_TARGET_DROP_BAR = 0.07

SIM_DOME_Z_THRESHOLD_MM = 591

# drip cut then dome cut
DOME_Z_THRESHOLD_CM_PRE = 6
DOME_Z_THRESHOLD_CM_POST = -3

# an exposure whose PT2121 live time is at or under this is dropped: a bubble that nucleates in the first
# moment of an exposure is not telling us about the threshold.
#
# Because the cut is on the live time ITSELF, the window in which a bubble could have been counted starts
# at the cut, not at zero -- an exposure that nucleated at 0.4 s is gone, so the first second of every
# surviving exposure is time in which nothing countable could have happened. Crediting it as exposure
# would inflate the denominator and pull every rate down. So a surviving exposure contributes
# (live time - LIVETIME_CUT_SEC), i.e. the rates count from the start of the cut.
#
# Applied identically to the source and background sides. It has to be: the two rates are subtracted, so
# if one counted from 0 and the other from 1 s they would be measuring different windows and the
# difference would carry that mismatch as a bias.
LIVETIME_CUT_SEC = 1.0


# alternate grouping: instead of one group per (p, T) setpoint, pool every (p, T) whose Seitz threshold
# falls inside the same alternateSeitzBinEdges range into a single group (pooled events, pooled livetime,
# pooled background, livetime-weighted Seitz threshold + sim prediction). Anything above the top edge is
# left alone, one group per (p, T) as usual. Writes to its own PLOTS_ROOT so it never overwrites the
# normal per-(p,T) output
useAlternateSeitzBinning = True
alternateSeitzBinEdges = [0.4, 0.6, 0.8, 1.0, 1.2, 1.4]

# granularity of the rebinned plots/fits: normally the 5 raw multiplicity classes collapse to 3 bins
# (1, 2, 3+); with this True they collapse to 2 instead (1, 2+), pooling every multi-bubble event into a
# single bin. Only touches the rebinned path -- useRebinnedThresholdPlots = False still draws all 5 raw
# classes either way. Writes to its own PLOTS_ROOT so it never overwrites the 3-bin output
useTwoPlusRebinning = False
rebinnedBinLabels = ["1", "2+"] if useTwoPlusRebinning else ["1", "2", "3+"]

# which multiplicity bins pin the shared sim-to-data normalization -- the source-activity number that puts
# the simulation on the same absolute scale as the data. "total" sums every bin, the historical behaviour.
# "singles" pins it on the multiplicity-1 bin alone: the highest-statistics bin, and the one least
# sensitive to how the sim models multi-bubble production, so it's the cleanest handle on the source rate
# itself -- at the cost of no longer forcing the summed rate to match. Either way the resulting number and
# what it does to every other bin is printed by report_source_activity_normalization()
# ABSOLUTE (source-activity) normalization -- the alternative to fitting the sim to the data at all.
# The sim generates a fixed number of neutrons into 4pi; the source emits a known number per minute; the
# ratio is how many minutes of running the simulation represents. Dividing raw sim counts by that gives a
# count/min rate with no reference to the measurement, so sim and data can be compared on absolute footing
# and any disagreement is a real disagreement rather than something the normalization has absorbed.
#
# Activity and neutron rate are Table I of the 2017 calibration paper, referenced to its stated date.
# Both alpha decay and spontaneous fission remove 252Cf nuclei, so the neutron rate (proportional to the
# atom count) decays with the TOTAL half-life, not the SF branch. Sanity check on the two table numbers:
# 0.0113 MBq = 11300 decays/s, x 3.092% SF branch x 3.7676 n/fission = 1316 n/s, vs the quoted 1320.
CF252_HALF_LIFE_YEARS = 2.645
CF252_REFERENCE_DATE = date(2017, 10, 11)
CF252_NEUTRONS_PER_SECOND_AT_REFERENCE = 1320.0
# when to evaluate the decayed activity. Sensitivity is roughly 2% per month, which dominates every other
# input here (half-life uncertainty and calendar convention are both sub-1%), so this is the number to
# revisit first if the absolute scale looks off -- ideally the midpoint of the neutron run period
ACTIVITY_EVALUATION_DATE = date(2026, 1, 1)
SIMULATED_NEUTRONS = 50e6


def cf252_neutrons_per_minute(evaluationDate=ACTIVITY_EVALUATION_DATE):
    years = (evaluationDate - CF252_REFERENCE_DATE).days / 365.25
    return CF252_NEUTRONS_PER_SECOND_AT_REFERENCE * 2 ** (-years / CF252_HALF_LIFE_YEARS) * 60


# minutes of source running the simulation represents, and the count/min-per-sim-count factor that follows.
# Computed here rather than further down because build_groups() needs it -- in this variant it IS the
# normalization, so there is nothing left to fit afterwards
simulatedMinutes = SIMULATED_NEUTRONS / cf252_neutrons_per_minute()
activityNormalizationFactor = 1.0 / simulatedMinutes

# "plotsFittedNorm" for a normal run, so this variant's output never overwrites neutronExcluded.py's own
# "plots/" -- see the module header. The two suffixes are kept because useAlternateSeitzBinning and
# useTwoPlusRebinning produce genuinely different output that would otherwise overwrite it
PLOTS_ROOT = (("plotsFittedNormAltBinning" if useAlternateSeitzBinning else "plotsFittedNorm")
              + ("Mult2Plus" if useTwoPlusRebinning else ""))

# rate tables land here rather than next to the script. Figures and the pipeline's own reports
# (summary.txt, perGroup.txt, chi2Debug.txt) stay under PLOTS_ROOT
OUTPUT_ROOT = os.path.join(SCRIPT_DIR, "output")

# every text file this script READS but never writes: the gray rate table write_final_rates() takes its
# alternate exposures from, and the per-run *_bubble_regions.txt files report_bubble_region_fractions()
# compares against. Nothing in here is ever modified
GRAY_INPUT_ROOT = os.path.join(SCRIPT_DIR, "gray-comparison-files")

# what write_final_rates() names the cf252-coffin-b_final_rates.txt-format tables it writes into
# OUTPUT_ROOT. The gray reference file they are modelled on is an INPUT, read from GRAY_INPUT_ROOT, and is
# named here partly so the write refuses to clobber it -- this script never edits that file
FINAL_RATES_FILENAME = "cf252-coffin-b-olivia.txt"
# the same table with the exposures taken from the reference file instead of this run's measured
# livetime -- see load_reference_exposures()
FINAL_RATES_TAU_FILENAME = "cf252-coffin-b-olivia-tauFromFile.txt"
FINAL_RATES_REFERENCE_FILENAME = "cf252-coffin-b_final_rates.txt"


def output_path(outputDir, filename):
    path = os.path.join(PLOTS_ROOT, outputDir, filename)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path


# a path inside OUTPUT_ROOT, creating it on first use
def final_output_path(filename):
    os.makedirs(OUTPUT_ROOT, exist_ok=True)
    return os.path.join(OUTPUT_ROOT, filename)


# Nothing in this pipeline writes to stdout -- every diagnostic, warning and caption number lands in one of
# three files under comparison/, so there are a few places to look instead of seven:
#
#   summary.txt   what the run produced: warnings, the sim-to-data normalization, how well sim and data
#                 line up, and each figure's fit statistics. The "read this first" file.
#   perGroup.txt  the per-(p,T)-group tables: setpoint pressures and the thresholds they imply, pre/post-cut
#                 rates, and any bins with no source-on events behind them.
#   chi2Debug.txt the verbose per-bin chi2 dump. Left on its own -- it is an order of magnitude longer than
#                 everything else combined and would bury the rest.
#
# Files open lazily on their first line, truncate once per run, and close at exit. Sections within a file
# appear in the order the pipeline produces them, each under a banner emitted once by report_section().
SUMMARY_REPORT = "summary.txt"
GROUP_REPORT = "perGroup.txt"

_openReports = {}
_reportSections = set()

def report_line(filename, line, outputDir="comparison"):
    handle = _openReports.get((outputDir, filename))
    if handle is None:
        handle = open(output_path(outputDir, filename), "w")
        _openReports[(outputDir, filename)] = handle
        atexit.register(handle.close)
    handle.write(line + "\n")


# append `lines` under `title` in a shared report file, emitting the banner only on the first call for that
# (file, title). Callers that build their whole body at once pass it in one go; callers that dribble a line
# at a time (caption stats, one per figure) can call repeatedly and land under a single banner
def report_section(filename, title, lines=(), outputDir="comparison"):
    key = (outputDir, filename, title)
    if key not in _reportSections:
        if (outputDir, filename) in _openReports:   # not the first section in this file
            report_line(filename, "", outputDir)
        _reportSections.add(key)
        report_line(filename, "=" * 78, outputDir)
        report_line(filename, title, outputDir)
        report_line(filename, "=" * 78, outputDir)
    for line in lines:
        report_line(filename, line, outputDir)


# every chi_squared_diagnostic(debug=True) line
def chi2_debug(line):
    report_line("chi2Debug.txt", line)


# excluded if the handscanner labeled it dome, or its precomputed position flag says it's past whichever dome-cut threshold is in play
def make_is_region_excluded(positionDomeFlags):
    def is_region_excluded(i, region):
        return SOURCE_LABELS.get(region) == "dome" or positionDomeFlags[i]
    return is_region_excluded


# True if a real reco hit's Z (detector frame) is past the dome cut
def is_dome_event(z_mm):
    return z_mm > DOME_Z_THRESHOLD_CM * 10


# same as is_dome_event(), but against an arbitrary threshold_cm instead of the module-level one
def is_dome_event_at(z_mm, threshold_cm):
    return z_mm > threshold_cm * 10


## combined reco + handscan region assignment, used by the FRACTION REPORTS ONLY
##
## The handscan gives one scan_source label per event, assigned by eye. The reconstruction gives a 3D
## position. Either one calling an event dome makes it dome; either one calling it wall makes it wall.
## Only single-bubble events get a region at all: the handscan labels a multi-bubble event as a whole,
## and multi-bubble events have no reliable 3D coord, so neither input means anything for them.
##
## Geometry is read off draw_detector_r2z_guides(), which draws the jar in r^2 vs z:
##   the cylindrical wall runs from z = -222.250 mm to z = -16.206 mm, inner surface r = 114.935 mm,
##   outer surface r = 120.015 mm (5.08 mm of glass), and the dome is the shoulder arc above that,
##   continuous with both wall surfaces at z = -16.206.
##
## "wall" is r >= WALL_RADIUS_MM, which also takes in anything reconstructed inside the glass itself,
## since that is further out still. "dome" reuses the pipeline's existing is_dome_event() threshold
## rather than a distance to the shoulder arc.
WALL_INNER_RADIUS_MM = 25.4 * 4.525
WALL_OUTER_RADIUS_MM = 25.4 * 4.725
# everything at or beyond this radius counts as wall
WALL_RADIUS_MM = 110.0
# how far inside the inner surface that puts the threshold, for the report text below
WALL_PROXIMITY_MM = WALL_INNER_RADIUS_MM - WALL_RADIUS_MM
CYLINDER_BOTTOM_MM = 25.4 * -8.75
CYLINDER_TOP_MM = 25.4 * (14.71997 - 15.358)

BULK_REGION_CODE = 0
WALL_REGION_CODE = 1
DOME_REGION_CODE = 2


# True if a reco coord sits at or beyond WALL_RADIUS_MM. Bounded in z by the cylinder itself: above
# CYLINDER_TOP_MM the boundary is the dome shoulder, not the wall, and that territory belongs to the
# dome test below
def is_reco_wall_event(coord):
    if coord is None:
        return False
    radius = float(np.hypot(coord[0], coord[1]))
    return radius >= WALL_RADIUS_MM and CYLINDER_BOTTOM_MM <= coord[2] <= CYLINDER_TOP_MM


# the single region code for one single-bubble event, OR-ing the handscan label with the reco position.
# An event can satisfy both tests -- the dome threshold sits 13.8 mm below the cylinder top, so the top
# of the cylindrical wall is inside it. Dome takes precedence there so the four regions stay mutually
# exclusive and still sum to the sample; report_*_fractions() counts the overlap separately
def combined_region(handscanRegion, coord):
    if handscanRegion == DOME_REGION_CODE or (coord is not None and is_dome_event(coord[2])):
        return DOME_REGION_CODE
    if handscanRegion == WALL_REGION_CODE or is_reco_wall_event(coord):
        return WALL_REGION_CODE
    return handscanRegion


# (isWall, isDome) before the precedence rule, so the overlap can be counted
def region_flags(handscanRegion, coord):
    isDome = handscanRegion == DOME_REGION_CODE or (coord is not None and is_dome_event(coord[2]))
    isWall = handscanRegion == WALL_REGION_CODE or is_reco_wall_event(coord)
    return isWall, isDome


# True if a sim hit's raw Z/mm is past the dome cut
def is_sim_dome_event(z_mm):
    return z_mm > SIM_DOME_Z_THRESHOLD_MM


## sim side

# same as cfconfBThresholds.get_multiplicity_counts(), but drops dome-region single-bubble events first when "dome" is in excludedRegionsOverride (defaults to module-level excludedRegions)
def get_multiplicity_counts(energy_threshold_kev, multiplicity_cut=MULTIPLICITY_CUT, excludedRegionsOverride=None):
    activeExcludedRegions = excludedRegions if excludedRegionsOverride is None else excludedRegionsOverride
    if "dome" not in activeExcludedRegions:
        return cfconf_get_multiplicity_counts(energy_threshold_kev, multiplicity_cut)

    df = get_event_coordinates(energy_threshold_kev)
    origMultiplicity = df.groupby("Event")["Event"].transform("size")
    isSingleBubble = origMultiplicity == 1
    isDomeHit = df["Z/mm"].apply(is_sim_dome_event)
    df = df.loc[~(isSingleBubble & isDomeHit), :]

    multiplicity = df.groupby("Event").size().tolist()
    if not multiplicity:
        # dome filter cut everything -- report zero counts instead of crashing on empty min()/max()
        zeros = np.zeros(multiplicity_cut)
        return zeros, zeros
    multiplicity_min = min(multiplicity)
    multiplicity_max = max(multiplicity)
    bin_num = multiplicity_max - multiplicity_min + 1
    bin_range = (multiplicity_min, multiplicity_max + 1)
    multiplicity_counts, _ = np.histogram(multiplicity, bins=bin_num, range=bin_range)

    ratio_sigma = (
        np.sqrt(multiplicity_counts) * multiplicity_counts[0]
        + multiplicity_counts * np.sqrt(multiplicity_counts[0])
    ) / multiplicity_counts[0] ** 2

    return multiplicity_counts[:multiplicity_cut], ratio_sigma[:multiplicity_cut]


## real-data side

SINGLE_BUBBLE_MULT = 1
MAX_FRAME = 50

# background rate calculation for subtraction
## warm annular
#backgroundRunsWarm = ["20251113_9","20251113_10","20251113_11","20251114_0","20251114_1","20251114_6","20251114_36","20251114_37","20251115_0","20251115_1","20251115_2","20251115_3","20251115_4","20251115_5","20251116_1","20251116_2","20251117_0","20251117_1","20251126_7","20251126_8","20251127_0","20251127_1","20251127_2","20251127_3","20251127_4","20251127_5","20251128_0","20251128_1","20251128_2","20251128_3","20251128_4","20251129_0","20251129_1","20251129_2","20251129_3","20251129_4","20251129_5","20251130_0","20251130_1","20251130_2","20251130_3","20251130_4","20251130_5",]
backgroundRunsWarm = []
## cold annular
backgroundRunsCold = ["20260117_0","20260117_1","20260117_2","20260117_3","20260117_4","20260118_0","20260118_1","20260118_2","20260118_3","20260118_4","20260119_0","20260119_1","20260119_2","20260120_0","20260120_1",]
## 199k
backgroundRunsHot = ["20260205_0","20260205_1","20260205_2","20260205_3","20260205_4","20260217_7","20260217_8","20260217_9","20260217_10","20260217_11","20260217_12","20260217_13","20260218_0","20260218_1","20260218_2","20260218_3","20260218_4","20260218_5","20260218_6","20260218_15","20260218_16","20260219_0","20260219_1","20260219_2","20260219_3","20260219_4","20260219_5","20260219_6","20260219_7","20260219_8","20260219_9","20260219_10","20260219_11","20260220_1","20260220_2","20260220_3","20260220_4",]
# ones to use for rate calculation
backgroundList = backgroundRunsWarm + backgroundRunsCold + backgroundRunsHot

# neutron source runs
# config A
## warm annular
neutronRunsWarm = ["20260107_3", "20260107_4", "20260107_5", "20260107_6", "20260107_7", "20260108_0", "20260108_1", "20260108_2", "20260108_3"]
## cold annular
neutronRunsCold = []
## 119K
neutronRunsHot = []

# config B
## warm annular
#neutronRunsWarmB = ["20260108_4", "20260108_5", "20260108_6", "20260108_7", "20260108_8", "20260109_0"]
neutronRunsWarmB = []
## cold annular
neutronRunsColdB = ["20260122_3","20260122_4","20260122_5","20260122_6","20260123_0","20260123_1","20260123_2","20260123_3","20260123_4","20260123_8","20260123_9","20260123_10","20260124_0","20260124_1","20260124_3","20260124_4","20260124_5","20260125_0","20260125_1","20260125_2","20260125_3","20260125_4","20260125_5","20260125_6","20260125_7","20260125_8"]
## 119K
neutronRunsHotB = ["20260205_12","20260205_13","20260205_14","20260205_15","20260205_16","20260205_17","20260205_18","20260206_0","20260206_1","20260206_2","20260206_3","20260206_4","20260206_5","20260206_6","20260206_7","20260213_1","20260213_2","20260213_3","20260213_4","20260213_5","20260213_6","20260213_7","20260213_8","20260213_9","20260214_0","20260214_1","20260214_2","20260214_3","20260214_4","20260214_5","20260214_6","20260214_7","20260214_8","20260214_9","20260214_10","20260214_11","20260214_12","20260214_13","20260214_14","20260215_0","20260215_1","20260215_2","20260215_3","20260215_4","20260215_5","20260215_6","20260215_7","20260215_8","20260215_9","20260215_10","20260215_11","20260215_12","20260215_13","20260215_14","20260216_0","20260216_1","20260216_2","20260216_3","20260216_4","20260216_5","20260216_6","20260216_7","20260216_8","20260216_9","20260216_10","20260216_11","20260216_12","20260216_13","20260217_0","20260217_1","20260217_2","20260217_3","20260217_4","20260217_5","20260217_6"]

## ones that are used for this graph
useConfigB = True

neutronRuns = (neutronRunsWarmB + neutronRunsColdB + neutronRunsHotB) if useConfigB \
    else (neutronRunsWarm + neutronRunsCold + neutronRunsHot)

## 119.6K runs, source and background alike -- run strings are full "date_run" (e.g. "20260205_12"), same form as the run column in the handscan files
hotRuns = neutronRunsHot + neutronRunsHotB + backgroundRunsHot

# nominal temperature in K of the run a given event came from
def run_temperature(run):
    return 119.6 if run in hotRuns else 116.7

## real-data region exclusion (scan_source codes); "dome" also gates the sim Z cut above
# scan_source code -> label, same mapping as combine_handscans.py / EventDisplay
SOURCE_LABELS = {0: "bulk", 1: "wall", 2: "dome", 3: "bellows", 4: "other"}
#excludedRegions = []
excludedRegions = ["dome"]

# returns (run, ev, mult, region) once per unique event whose run is in runList
def iter_matched_events(dirpath, runList):
    checked = set()
    for path in glob.glob(os.path.join(dirpath, "*.txt")):
        with open(path, 'r', encoding='utf-8') as f:
            for raw in f:
                parts = raw.split()
                if len(parts) < 5:
                    continue
                run, ev = parts[0], parts[1]
                try:
                    mult = int(parts[4])
                    region = int(parts[3])
                except ValueError:
                    continue
                if (ev, run) in checked:
                    continue
                checked.add((ev, run))
                if run in runList:
                    yield run, ev, mult, region


# confirmed single-bubble (mult == 1) events, for the reco-position lookup below
def iter_single_bubble_events(dirpath, runList):
    for run, ev, mult, region in iter_matched_events(dirpath, runList):
        if mult == SINGLE_BUBBLE_MULT:
            yield run, ev, region


def _is_bad_coord(coord):
    return np.isnan(coord).any() or coord[0] <= -999


# {(ev, frame): first valid 3D coord}, handles both reco.sbc "frame" shapes seen so far: (N,50) blocks and plain 1D per-row values
def build_reco_lookup(reconInfo):
    recoLookup = {}
    reconEv = reconInfo["ev"]
    reconFrame = reconInfo["frame"]
    coords3D = reconInfo["coords_3D"]
    nCoords = len(coords3D)
    for i in range(len(reconEv)):
        ev_i = reconEv[i]
        frameRow = reconFrame[i]
        if np.ndim(frameRow) == 0:
            frameNumIdxPairs = [(frameRow, i)]
        else:
            frameNumIdxPairs = [(frameRow[j], i + j) for j in range(len(frameRow))]
        for frameNum, idx in frameNumIdxPairs:
            if idx >= nCoords:
                continue
            key = (ev_i, frameNum)
            if key in recoLookup:
                continue
            coord = coords3D[idx]
            if not _is_bad_coord(coord):
                recoLookup[key] = coord
    return recoLookup


# first non-error 3D coord for evNum, scanning frames in ascending order
def first_valid_coord(recoLookup, evNum):
    evNum = int(evNum)
    for f in range(MAX_FRAME):
        coord = recoLookup.get((evNum, f))
        if coord is not None:
            return coord
    return None


# excluded if the handscanner labeled it dome (region == 2) or its reco Z is past the dome cut
def is_real_dome_event(region, coord):
    if region == 2:
        return True
    return coord is not None and is_dome_event(coord[2])


# (run, ev, x, y, z) per confirmed single-bubble event; drops dome events unless applyDomeExclusion=False (then kept, only counted for the print below). label is just for that print
def load_single_bubble_positions(runList, label="real data", applyDomeExclusion=True):
    positions = []
    recoLookupCache = {}
    regionExcludedCount = 0
    zExcludedCount = 0
    for run, ev, region in iter_single_bubble_events(HANDSCAN_DIR, runList):
        if run not in recoLookupCache:
            recoPath = os.path.join(RECON_DIR, run, "reco.sbc")
            if not os.path.exists(recoPath):
                recoLookupCache[run] = None
            else:
                recoLookupCache[run] = build_reco_lookup(Streamer(recoPath).to_dict())
        recoLookup = recoLookupCache[run]
        if recoLookup is None:
            continue
        coord = first_valid_coord(recoLookup, ev)
        if coord is None:
            continue
        if "dome" in excludedRegions and is_real_dome_event(region, coord):
            if region == 2:
                regionExcludedCount += 1
            else:
                zExcludedCount += 1
            if applyDomeExclusion:
                continue
        positions.append((run, ev, float(coord[0]), float(coord[1]), float(coord[2])))
    return positions


# (pset_lo, pset_hi, temp) for one event, read from its run's event.sbc
def _event_pset_temp(run, ev, evDataCache):
    if run not in evDataCache:
        evPath = os.path.join(RECON_DIR, run, "event.sbc")
        evDataCache[run] = Streamer(evPath).to_dict() if os.path.exists(evPath) else None
    evData = evDataCache[run]
    if evData is None:
        return None
    for i in range(len(evData["ev"])):
        if int(evData["ev"][i]) == int(ev):
            return float(evData["pset_lo"][i]), float(evData["pset_hi"][i]), run_temperature(run)
    return None


# background bin counts (multiplicity 1,2,3,4,5+) and live time, split by the (pressure, temperature) each event was taken at,
# so every source group is subtracted with background from its own operating point. Returns {(p, T): (binCounts, liveTime)};
# the live time one exposure contributes, measured from the start of the cut, or None if the cut drops the
# exposure entirely. The single place LIVETIME_CUT_SEC is applied -- see it for why the offset is there
def cut_livetime(rawLiveTime):
    liveTime = float(rawLiveTime)
    return liveTime - LIVETIME_CUT_SEC if liveTime > LIVETIME_CUT_SEC else None


# background events without a fixed setpoint belong to no group and are dropped, same rule pToUse applies to the source side
def load_background():
    subdirs = {os.path.basename(p.rstrip(os.sep)) for p in glob.glob(os.path.join(RECON_DIR, '*/'))}
    byPT = {}
    # same events, same livetime, binned by combined_region() (handscan label OR reco position -- see
    # final_rates_source_by_pt()) instead of multiplicity -- only write_final_rates() wants this, so it
    # rides along here rather than costing a second pass over every exposure.sbc/event.sbc. Only
    # single-bubble events get a reco coord (below), so a multi-bubble event's region is its raw handscan
    # label unchanged, same as the source side. {(p, T): [count per SOURCE_LABELS code]}, livetime lives in
    # byPT. This does its own reco.sbc lookups rather than reusing background_combined_regions()'s --
    # that function is single-bubble-only and un-livetime-cut, built for the (different) fraction reports,
    # so its counts are not interchangeable with the livetime-cut, all-multiplicity ones this needs
    regionByPT = {}
    # the same thing restricted to single-bubble events. Needed separately because gray's N Localized is
    # single-bubble by construction, so the cross-comparison in wallOrDomeFractions.txt cannot use the
    # all-multiplicity counts above without comparing two different samples
    singleRegionByPT = {}
    # parallel per-event arrays, the same shape load_neutron_events() returns for the source side
    # (bubbleCount/sourceTimes/psetsTemps/runEvs). Built so the background side can go through the exact
    # same bin_multiplicities()/compute_position_dome_flags() machinery the source side uses for its dome
    # cut, rather than a hand-rolled equivalent that could drift out of sync with it. Every matched event
    # with a resolvable pset goes in, stable or not -- same as the source side, which lets inGroup() do the
    # stability filtering downstream instead of filtering here
    bubbleCountEvents, backgroundTimes, psetsTempsEvents, runEvs = [], [], [], []
    expDataCache, evDataCache, recoLookupCache = {}, {}, {}
    for run, ev, mult, region in iter_matched_events(HANDSCAN_DIR, backgroundList):
        if run not in subdirs:
            continue
        if run not in expDataCache:
            expDataCache[run] = Streamer(os.path.join(RECON_DIR, run, 'exposure.sbc')).to_dict()
        expData = expDataCache[run]
        for i in range(len(expData["ev"])):
            if int(expData["ev"][i]) == int(ev):
                cutLiveTime = cut_livetime(expData['PT2121_livetime'][i])
                if cutLiveTime is None:
                    break
                pset = _event_pset_temp(run, ev, evDataCache)
                if pset is None:
                    break
                bubbleCountEvents.append((mult, region))
                runEvs.append((run, ev))
                backgroundTimes.append(cutLiveTime)
                psetsTempsEvents.append(pset)
                if pset[0] == pset[1]:
                    binCounts, liveTime = byPT.setdefault((pset[0], pset[2]), ([0] * 5, 0.0))
                    regionCounts = regionByPT.setdefault((pset[0], pset[2]), [0] * len(SOURCE_LABELS))
                    singleCounts = singleRegionByPT.setdefault((pset[0], pset[2]),
                                                               [0] * len(SOURCE_LABELS))
                    if mult != 0:
                        binCounts[min(mult, 5) - 1] += 1
                        coord = None
                        if mult == SINGLE_BUBBLE_MULT:
                            if run not in recoLookupCache:
                                recoPath = os.path.join(RECON_DIR, run, "reco.sbc")
                                recoLookupCache[run] = build_reco_lookup(Streamer(recoPath).to_dict()) \
                                    if os.path.exists(recoPath) else None
                            recoLookup = recoLookupCache[run]
                            if recoLookup is not None:
                                coord = first_valid_coord(recoLookup, ev)
                        combinedRegion = combined_region(region, coord)
                        if 0 <= combinedRegion < len(regionCounts):
                            regionCounts[combinedRegion] += 1
                            if mult == SINGLE_BUBBLE_MULT:
                                singleCounts[combinedRegion] += 1
                    byPT[(pset[0], pset[2])] = (binCounts, liveTime + cutLiveTime)
                break
    return byPT, regionByPT, singleRegionByPT, bubbleCountEvents, backgroundTimes, psetsTempsEvents, runEvs

# raw PT1101 transducer reading -> LAr pressure at the transducer [bar]
def corrected_pt1101(rawPressure):
    return float(rawPressure) * PT1101_GAIN - PT1101_OFFSET


# ... and down to the LAr pressure in the target volume, which is what the Seitz model wants. NaN in, NaN
# out: a run whose slow DAQ didn't record PT1101 is dropped in build_pressure_distributions() rather than
# poisoning that setpoint's mean
def corrected_pt1101_target(rawPressure):
    return corrected_pt1101(rawPressure) - PT1101_TARGET_DROP_BAR


# (mult, region) per neutron-run event, its live time, its target-volume PT1101 pressure, and its (pset_lo, pset_hi, temp)
def load_neutron_events():
    bubbleCount, sourceTimes, sourcePressures, psetsTemps, runEvs = [], [], [], [], []
    expDataCache, evDataCache = {}, {}
    for run, ev, mult, region in iter_matched_events(HANDSCAN_DIR, neutronRuns):
        if run not in expDataCache:
            expDataCache[run] = Streamer(f'{RECON_DIR}{run}/exposure.sbc').to_dict()
        expData = expDataCache[run]
        liveTime, pressure = None, float("nan")
        for i in range(len(expData["ev"])):
            if int(expData["ev"][i]) == int(ev):
                # None here is the LIVETIME_CUT_SEC cut firing, and it drops the exposure below exactly the
                # way a missing one does. This side used to take the raw live time with no cut at all,
                # which left it counting from 0 while the background counted from 1 s
                liveTime = cut_livetime(expData['PT2121_livetime'][i])
                if 'PT1101_pressure' in expData:
                    pressure = corrected_pt1101_target(expData['PT1101_pressure'][i])
                break
        pset = _event_pset_temp(run, ev, evDataCache)
        # all five lists are indexed in parallel downstream, so an event the exposure or event file doesn't
        # have -- or one the live-time cut drops -- goes in none of them: appending to only some would
        # shift every later event's livetime/(p, T)
        if liveTime is None or pset is None:
            continue
        bubbleCount.append((mult, region))
        runEvs.append((run, ev))
        sourceTimes.append(liveTime)
        sourcePressures.append(pressure)
        psetsTemps.append(pset)
    return bubbleCount, sourceTimes, sourcePressures, psetsTemps, runEvs


# per neutron-run event, whether it's a confirmed single-bubble event past threshold_cm -- multi-bubble events get no reliable 3D coord, so they're always False here
def compute_position_dome_flags(bubbleCount, runEvs, threshold_cm):
    flags = [False] * len(bubbleCount)
    recoLookupCache = {}
    for i, (mult, region) in enumerate(bubbleCount):
        if mult != SINGLE_BUBBLE_MULT:
            continue
        run, ev = runEvs[i]
        if run not in recoLookupCache:
            recoPath = os.path.join(RECON_DIR, run, "reco.sbc")
            recoLookupCache[run] = build_reco_lookup(Streamer(recoPath).to_dict()) if os.path.exists(recoPath) else None
        recoLookup = recoLookupCache[run]
        if recoLookup is None:
            continue
        coord = first_valid_coord(recoLookup, ev)
        if coord is None:
            continue
        flags[i] = is_dome_event_at(coord[2], threshold_cm)
    return flags

# reco coord (or None) per neutron-run event, parallel to bubbleCount. Only single-bubble events get one:
# a multi-bubble event has no reliable 3D coord, and the fraction reports only assign regions to single
# bubbles anyway. Same reco lookup compute_position_dome_flags() uses, walked once
def compute_single_bubble_coords(bubbleCount, runEvs):
    coords = [None] * len(bubbleCount)
    recoLookupCache = {}
    for i, (mult, _region) in enumerate(bubbleCount):
        if mult != SINGLE_BUBBLE_MULT:
            continue
        run, ev = runEvs[i]
        if run not in recoLookupCache:
            recoPath = os.path.join(RECON_DIR, run, "reco.sbc")
            recoLookupCache[run] = build_reco_lookup(Streamer(recoPath).to_dict()) \
                if os.path.exists(recoPath) else None
        recoLookup = recoLookupCache[run]
        if recoLookup is not None:
            coords[i] = first_valid_coord(recoLookup, ev)
    return coords


# counts per multiplicity bin (1..5+) and total live time. The two predicates do different jobs, and
# conflating them is what the old single `keep` got wrong:
#
#   inGroup(i)            is this exposure part of this group at all? A (p, T) outside the group never
#                         ran for it, so it contributes neither live time nor count.
#   isVetoed(i, region)   should this exposure's BUBBLE be thrown out (dome region, position cut)?
#
# A vetoed exposure keeps its live time. The detector really was live for it -- we are declining to count
# the bubble it produced, not claiming the exposure never happened. Dropping the live time too would
# remove real exposure from the denominator and inflate every rate, and by exactly the wrong amount: the
# more dome events a setpoint has, the more exposure it would lose.
def bin_multiplicities(bubbleCount, sourceTimes, inGroup, isVetoed):
    binCounts = [0] * 5
    sourceTime = 0.0
    for i, (mult, region) in enumerate(bubbleCount):
        if not inGroup(i):
            continue
        sourceTime += sourceTimes[i]
        if mult != 0 and not isVetoed(i, region):
            binCounts[min(mult, 5) - 1] += 1
    return binCounts, sourceTime


# exact one-sided 68% Poisson upper limit for an observed count of 0: P(0 | mu) = e^-mu = 1 - 0.68
ZERO_COUNT_UPPER = -np.log(1 - 0.68)

# background-subtracted counts and their (asymmetric) errors, background scaled to sourceTime
def background_subtract(binCounts, sourceTime, backgroundBinCounts, backgroundTime):
    scale = sourceTime / backgroundTime if backgroundTime > 0 else 0.0
    backBins = [b * scale for b in backgroundBinCounts]
    # an empty Poisson bin can't have fluctuated down, but it can have fluctuated up from a true mean as
    # large as -ln(1 - 0.68) = 1.14 (the exact one-sided 68% upper limit for N = 0), so zero-count bins get
    # errLow = 0 and that upper limit instead of sqrt(0) = 0 on both sides. Applied to the source counts and
    # the background counts alike -- doing it for only one of them is what left an empty source bin with no
    # upper error at all while still carrying a (downward, unphysical) background-driven lower error
    backErrorLow = [np.sqrt(b) * scale if b > 0 else 0.0 for b in backgroundBinCounts]
    backErrorHigh = [(np.sqrt(b) if b > 0 else ZERO_COUNT_UPPER) * scale for b in backgroundBinCounts]
    binCountErrorLow = [np.sqrt(c) if c > 0 else 0.0 for c in binCounts]
    binCountErrorHigh = [np.sqrt(c) if c > 0 else ZERO_COUNT_UPPER for c in binCounts]
    backSubBins = [c - b for c, b in zip(binCounts, backBins)]
    # background high -> subtracted rate pulled down; background low -> subtracted rate pulled up
    backSubErrorLow = [np.sqrt(countErr**2 + backErr**2) for countErr, backErr in zip(binCountErrorLow, backErrorHigh)]
    backSubErrorHigh = [np.sqrt(countErr**2 + backErr**2) for countErr, backErr in zip(binCountErrorHigh, backErrorLow)]
    return (backBins, backErrorLow, backErrorHigh, backSubBins, backSubErrorLow, backSubErrorHigh,
            binCountErrorLow, binCountErrorHigh)


# PLOTTING ONLY -- a negative rate is unphysical, so the drawn point sits at 0 with no lower error and an upper error
# trimmed to the 68% upper edge the fluctuation actually allows. Deliberately NOT fed to the fits: flooring an observed
# value pulls it toward the prediction, and the trimmed error becomes the chi2 denominator, so a deficit bin either
# blows up or (once its upper edge reaches 0) drops out of the sum entirely. Statistics use the unfloored values below
def floor_at_zero(values, errLow, errHigh):
    flooredErrHigh = [max(v + e, 0.0) if v < 0 else e for v, e in zip(values, errHigh)]
    flooredErrLow = [min(e, v) if v > 0 else 0.0 for v, e in zip(values, errLow)]
    return [max(v, 0.0) for v in values], flooredErrLow, flooredErrHigh


# collapse the 5 multiplicity classes (1,2,3,4,5+) down to rebinnedBinLabels: [1, 2, 3+], or [1, 2+] with
# useTwoPlusRebinning. Every class past the last label is summed into it
def rebin(values):
    keptBins = len(rebinnedBinLabels) - 1
    return list(values[:keptBins]) + [sum(values[keptBins:])]


def counts_to_ratios(counts):
    total = sum(counts)
    return [c / total for c in counts]

# scale a multiplicity ratio shape to match a target total
def seitz_count(ratios, total):
    scale = total / sum(ratios)
    return [scale * r for r in ratios]

# load in data (excludedRegions-independent -- region filtering happens per-group below)
(backgroundByPT, backgroundRegionsByPT, backgroundSingleRegionsByPT,
 backgroundBubbleCount, backgroundTimes, backgroundPsetsTemps, backgroundRunEvs) = load_background()


# {(p, T): [count per region]} for single-bubble BACKGROUND events under the combined reco+handscan
# assignment, plus the same diagnostics the source side reports. load_background() now keeps run/ev too
# (backgroundRunEvs), but this stays a separate pass: it is single-bubble-only and NOT livetime-cut, built
# for the (different) fraction reports, so its counts are not interchangeable with load_background()'s
# livetime-cut, all-multiplicity ones
def background_combined_regions():
    byPT = {}
    diagnostics = {"single": 0, "noCoord": 0, "bothFlagged": 0, "movedByReco": 0}
    recoLookupCache, evDataCache = {}, {}
    for run, ev, region in iter_single_bubble_events(HANDSCAN_DIR, backgroundList):
        pset = _event_pset_temp(run, ev, evDataCache)
        if pset is None or pset[0] != pset[1]:
            continue
        if run not in recoLookupCache:
            recoPath = os.path.join(RECON_DIR, run, "reco.sbc")
            recoLookupCache[run] = build_reco_lookup(Streamer(recoPath).to_dict()) \
                if os.path.exists(recoPath) else None
        recoLookup = recoLookupCache[run]
        coord = first_valid_coord(recoLookup, ev) if recoLookup is not None else None

        diagnostics["single"] += 1
        if coord is None:
            diagnostics["noCoord"] += 1
        isWall, isDome = region_flags(region, coord)
        if isWall and isDome:
            diagnostics["bothFlagged"] += 1
        combined = combined_region(region, coord)
        if combined != region:
            diagnostics["movedByReco"] += 1

        counts = byPT.setdefault((pset[0], pset[2]), [0] * len(SOURCE_LABELS))
        if 0 <= combined < len(counts):
            counts[combined] += 1
    return byPT, diagnostics


backgroundCombinedRegionsByPT, backgroundCombinedDiagnostics = background_combined_regions()


# the same for the SOURCE side, off the parallel arrays that already carry pset and reco coord
def source_combined_regions():
    byPT = {}
    diagnostics = {"single": 0, "noCoord": 0, "bothFlagged": 0, "movedByReco": 0}
    for (mult, region), pset, coord in zip(bubbleCount, psetsTemps, sourceSingleBubbleCoords):
        if mult != SINGLE_BUBBLE_MULT or float(pset[0]) != float(pset[1]):
            continue
        diagnostics["single"] += 1
        if coord is None:
            diagnostics["noCoord"] += 1
        isWall, isDome = region_flags(region, coord)
        if isWall and isDome:
            diagnostics["bothFlagged"] += 1
        combined = combined_region(region, coord)
        if combined != region:
            diagnostics["movedByReco"] += 1

        counts = byPT.setdefault((float(pset[0]), float(pset[2])), [0] * len(SOURCE_LABELS))
        if 0 <= combined < len(counts):
            counts[combined] += 1
    return byPT, diagnostics


# background counts/live time for one (p, T), or summed over several for the all-groups-averaged pipeline.
# isVetoed applies the SAME dome exclusion the source side's bin_multiplicities() call gets for this
# pipeline (bubble vetoed, live time kept) -- see the module header note on bin_multiplicities() for why a
# vetoed exposure keeps its live time. Routed through bin_multiplicities() itself, on the background's own
# parallel arrays (backgroundBubbleCount/backgroundTimes/backgroundPsetsTemps), rather than a hand-rolled
# equivalent, so the two sides cannot drift out of sync with each other
def background_for(pTs, isVetoed):
    memberSet = set(pTs)
    return bin_multiplicities(
        backgroundBubbleCount, backgroundTimes,
        inGroup=lambda i: (backgroundPsetsTemps[i][0], backgroundPsetsTemps[i][2]) in memberSet,
        isVetoed=isVetoed,
    )

bubbleCount, sourceTimes, sourcePressures, psetsTemps, neutronRunEvs = load_neutron_events()

# reco-Z-based dome flags for the withoutDomeCut/withDomeCut pipeline: "pre" = always-applied loose threshold, "post" = tighter threshold on top
# reco positions for the combined region assignment the fraction reports use
sourceSingleBubbleCoords = compute_single_bubble_coords(bubbleCount, neutronRunEvs)
# defined further up next to its background twin, but only callable now that the arrays it reads exist
sourceCombinedRegionsByPT, sourceCombinedDiagnostics = source_combined_regions()

positionDomeFlagsPre = compute_position_dome_flags(bubbleCount, neutronRunEvs, DOME_Z_THRESHOLD_CM_PRE)
positionDomeFlagsPost = compute_position_dome_flags(bubbleCount, neutronRunEvs, DOME_Z_THRESHOLD_CM_POST)

# the background sample's own position dome flags, same thresholds, so build_groups() can dome-cut the
# background the same way it dome-cuts the source for each pipeline instead of leaving it uncut
backgroundPositionDomeFlagsPre = compute_position_dome_flags(backgroundBubbleCount, backgroundRunEvs,
                                                              DOME_Z_THRESHOLD_CM_PRE)
backgroundPositionDomeFlagsPost = compute_position_dome_flags(backgroundBubbleCount, backgroundRunEvs,
                                                               DOME_Z_THRESHOLD_CM_POST)

# (p, T) pairs with a fixed pressure setpoint, deduped
pToUse = sorted({(float(lo), float(t)) for lo, hi, t in psetsTemps if float(lo) == float(hi)})


# every per-event target-volume PT1101 pressure [bar] recorded at each (pset, temperature), as a plain
# list -- the distribution whose median becomes that setpoint's Seitz pressure. Events with no PT1101
# reading are skipped; a setpoint with no valid reading at all is simply absent from the map
def build_pressure_distributions():
    byPT = {}
    for pressure, liveTime, pset in zip(sourcePressures, sourceTimes, psetsTemps):
        if not np.isfinite(pressure):
            continue
        byPT.setdefault((float(pset[0]), float(pset[2])), []).append((pressure, liveTime))
    return byPT


pressureDistributionByPT = build_pressure_distributions()
# median rather than mean: each setpoint's distribution is a tight spike at the operating pressure plus a
# sparse high tail of exposures that never really expanded (see report_pressure_distributions below), and
# that tail pulls the mean several tenths of a bar off the spike while leaving the median on it
measuredPressureByPT = {pT: float(np.median([pressure for pressure, _ in samples]))
                        for pT, samples in pressureDistributionByPT.items()}

# nominal-pset setpoints with no usable PT1101 reading -- those fall back to the pset itself below, so
# say so loudly rather than silently mixing measured and nominal pressures across groups
missingPressurePTs = [pT for pT in pToUse if pT not in measuredPressureByPT]
if missingPressurePTs:
    report_section(SUMMARY_REPORT, "Run warnings",
                   [f"WARNING: no valid PT1101 reading for {len(missingPressurePTs)} setpoint(s) "
                    f"{missingPressurePTs} -- falling back to the nominal pset for their Seitz threshold"])


# Seitz threshold [keV] for one (P [bar], T [K]) setpoint. The (p, T) pair identifies the setpoint, but the
# pressure actually fed to the Seitz model is the median of that setpoint's target-volume PT1101
# distribution (see build_pressure_distributions), not the nominal pset
def seitz_model(p, T):
    pressure = measuredPressureByPT.get((float(p), float(T)), p)
    return sm.SeitzModel(pressure * 14.5038, -273.15 + T, 'argon')


def seitz_threshold(p, T):
    return seitz_model(p, T).Q


# the distributions the medians above come from, printed and drawn one panel per setpoint. Worth
# eyeballing: the median only describes the operating point if the distribution really is one spike plus a
# tail, so a genuinely broad or bimodal setpoint (pressure ramps, a mislabelled pset) shows up here before
# it quietly shifts a group's threshold. Each panel marks the median fed to SeitzModel, the mean it was
# chosen over, and the nominal pset
def report_pressure_distributions(savepath):
    title = "Target-volume PT1101 pressure per setpoint (the median feeds SeitzModel)"
    lines = [f"    raw PT1101 * {PT1101_GAIN} - {PT1101_OFFSET} = pressure at the transducer, "
             f"minus {PT1101_TARGET_DROP_BAR} bar = target volume"]
    panels = [(pT, pressureDistributionByPT[pT]) for pT in pToUse if pT in pressureDistributionByPT]
    if not panels:
        lines.append("  none -- no setpoint has a usable PT1101 reading, every threshold falls back to its pset")
        report_section(GROUP_REPORT, title, lines)
        return

    for (pset, T), samples in panels:
        arr = np.asarray([pressure for pressure, _ in samples])
        expanded = np.asarray([pressure for pressure, liveTime in samples if liveTime > 1])
        nominalQ = sm.SeitzModel(pset * 14.5038, -273.15 + T, 'argon').Q
        lines.append(f"  pset = {pset:0.2f} bar, T = {T:0.1f} K:  n = {len(arr)} events, "
                     f"median = {np.median(arr):0.4f} bar, mean = {arr.mean():0.4f}, "
                     f"std = {arr.std(ddof=1) if len(arr) > 1 else 0.0:0.4f}, "
                     f"min = {arr.min():0.4f}, max = {arr.max():0.4f}")
        lines.append(f"      median - pset = {np.median(arr) - pset:+0.4f} bar  ->  "
                     f"Seitz {nominalQ:0.3f} keV (at pset) -> {seitz_threshold(pset, T):0.3f} keV (at median)")
        # exposures that never really expanded (the >1 s livetime cut load_background() already applies)
        # sit far above the setpoint. They're what the median is there to ignore, so show that it does:
        # a median that barely moves when they're dropped means the tail isn't setting the threshold
        if len(expanded) < len(arr):
            medianExpanded = np.median(expanded) if len(expanded) else float("nan")
            lines.append(f"      {len(arr) - len(expanded)} of {len(arr)} events have livetime <= 1 s; "
                         f"dropping them moves the median {np.median(arr):0.4f} -> {medianExpanded:0.4f} bar "
                         f"(the mean would move {arr.mean():0.4f} -> "
                         f"{expanded.mean() if len(expanded) else float('nan'):0.4f})")

    report_section(GROUP_REPORT, title, lines)

    nCols = min(3, len(panels))
    nRows = int(np.ceil(len(panels) / nCols))
    fig, axes = plt.subplots(nRows, nCols, figsize=(4.5 * nCols, 3.2 * nRows), squeeze=False)
    for ax, ((pset, T), samples) in zip(axes.ravel(), panels):
        arr = np.asarray([pressure for pressure, _ in samples])
        ax.hist(arr, bins=min(40, max(5, len(arr) // 5)), color="steelblue", edgecolor="white")
        ax.axvline(np.median(arr), color="red", linewidth=2,
                   label=f"median = {np.median(arr):0.3f} bar\n(Seitz {seitz_threshold(pset, T):0.2f} keV)")
        ax.axvline(arr.mean(), color="darkorange", linestyle=":", linewidth=2,
                   label=f"mean = {arr.mean():0.3f} bar (unused)")
        ax.axvline(pset, color="gray", linestyle="--", linewidth=1.5, label=f"pset = {pset:0.2f} bar")
        ax.set_title(f"pset {pset:0.2f} bar, {T:0.1f} K  (n = {len(arr)})", fontsize=10)
        ax.set_xlabel("target-volume PT1101 pressure [bar]", fontsize=9)
        ax.set_ylabel("events", fontsize=9)
        ax.tick_params(labelsize=8)
        ax.legend(fontsize=7)
    for ax in axes.ravel()[len(panels):]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(savepath)
    plt.close(fig)


report_pressure_distributions(output_path("comparison", "pt1101TargetPressureDistributions.png"))


## ---------------------------------------------------------------------------------------------------
## cf252-coffin-b_final_rates.txt-format output
##
## One row per (pset, temperature) in pToUse, in the same 26 columns and the same numpy.savetxt encoding as
## the reference file: a tab-joined "# "-prefixed header, then %.18e fields separated by single spaces.
##
## This table is deliberately NOT the multiplicity pipeline's numbers re-arranged. It is a straight
## event-count-over-exposure rate per setpoint, split by region rather than by bubble multiplicity, and
## with no dome CUT of any kind -- the Dome column only means anything if dome events are still in the
## sample. Nothing here feeds the plots or the fits.
##
## Written twice, from the same counts, differing only in the exposure the counts are divided by:
##   FINAL_RATES_FILENAME      this run's own measured PT2121 livetime, one per setpoint
##   FINAL_RATES_TAU_FILENAME  the exposure backed out of the reference file's own rate/error pairs
## See load_reference_exposures() for what the second one is and why it is per column rather than per row.
##
## NOTHING here is cut. The counts come straight off load_neutron_events(), so no dome cut (neither the
## handscan region == 2 exclusion nor the reco-Z position cut), no excludedRegions, no multiplicity cut --
## none of that machinery is on this path at all. Rate is every event recorded at the setpoint.
##
## Region mapping: combined_region() -- handscan label OR reco position, same rule the fraction reports
## use -- assigns each event a SOURCE_LABELS code, and codes 0/1/2/3 become the Bulk/Wall/Dome/Bottom
## columns (3 is the bellows, at the bottom of the vessel). It is NOT the raw handscan label alone: on a
## real run roughly half of single-bubble events get moved to a different region once reco position is
## consulted (report_bubble_region_fractions()'s diagnostics), and a raw-handscan-only split disagreed with
## the reference file's own (reco-based) region columns by an order of magnitude -- see
## report_wall_or_dome_fractions()'s cross-comparison table, which states outright that gray's columns are
## reco-based. Code 4 ("other") has no column in this format, but its events are still in Rate, because
## Rate is uncut -- so Rate >= the sum of the four region columns, strictly greater wherever the combined
## region is still "other" (or an unrecognized code) after the reco check -- unlike the reference file,
## whose Rate happens to equal its region sum exactly. write_final_rates() reports the difference per run
## so it is never a silent discrepancy.
FINAL_RATES_REGIONS = [("Bulk", 0), ("Wall", 1), ("Dome", 2), ("Bottom", 3)]
FINAL_RATES_UNCOLUMNED_REGIONS = [code for code in SOURCE_LABELS
                                  if code not in {code for _, code in FINAL_RATES_REGIONS}]
SECONDS_PER_HOUR = 3600.0

# columns the reference file carries between the Seitz threshold and the rates. Which SeitzModel attributes
# they come from is not pinned down anywhere in this repo and nothing here uses anything but .Q, so they are
# written as nan rather than guessed at -- the column is present and the format matches, the number is just
# absent. To fill one in, give it a "attribute name" second element and it gets read off the model object
# and written through unconverted, so check the units first: the reference file's headers say keV / nm /
# g/cc, and its numbers (rho_l ~ 1.18 g/cc for LAr near 117 K) say the raw attributes are already in those
FINAL_RATES_MODEL_FIELDS = [
    ("Eion [keV]", None),
    ("rion [nm]", None),
    ("rho_l [g/cc]", None),
]

# the five (Rate, Rate Error) pairs, in the order they appear in the format: the summed Rate first, then one
# per region. Index k here is the same k load_reference_exposures() and write_final_rates() use
FINAL_RATES_RATE_COLUMNS = ["Total"] + [name for name, _ in FINAL_RATES_REGIONS]
# where the first of those pairs starts: past Pressure/Temperature/Seitz Threshold and the model fields
FINAL_RATES_FIRST_RATE_COLUMN = 3 + len(FINAL_RATES_MODEL_FIELDS)


def _final_rates_model_field(model, attribute):
    if attribute is None:
        return float("nan")
    return float(getattr(model, attribute))


# {(p, T): ([count per SOURCE_LABELS code], total count, livetime [s])} for the source runs, uncut. The
# total is every bubble event at the setpoint whatever its region -- including "other" and including any
# scan_source code outside SOURCE_LABELS entirely -- so it is genuinely everything that was recorded, not
# the sum of the four columned regions.
#
# The livetime is every event's at that setpoint, region-independent -- the exposure is the same regardless
# of where in the vessel the bubble ended up, so all the rate columns are rates over that one livetime.
# mult == 0 means an exposure that produced no bubble: it contributes its livetime and no count.
#
# The region a count lands in is combined_region()'s (handscan label OR reco position, same as the fraction
# reports), not the raw handscan label alone -- sourceSingleBubbleCoords is None for every multi-bubble
# event, so those fall back to the raw label untouched; only single-bubble events can be moved by reco.
# This matters: on a real run roughly half of single-bubble events get moved by the reco check (see
# report_bubble_region_fractions()'s diagnostics), so a raw-handscan-only split reads nothing like gray's
# own region columns, which are reco-based -- see report_wall_or_dome_fractions()'s cross-comparison table.
def final_rates_source_by_pt():
    byPT = {}
    for (mult, region), liveTime, pset, coord in zip(bubbleCount, sourceTimes, psetsTemps,
                                                      sourceSingleBubbleCoords):
        if float(pset[0]) != float(pset[1]):
            continue
        key = (float(pset[0]), float(pset[2]))
        counts, total, liveTimeSec = byPT.setdefault(key, ([0] * len(SOURCE_LABELS), 0, 0.0))
        if mult != 0:
            total += 1
            combinedRegion = combined_region(region, coord)
            if 0 <= combinedRegion < len(counts):
                counts[combinedRegion] += 1
        byPT[key] = (counts, total, liveTimeSec + liveTime)
    return byPT


# raw count -> (rate [/hr], error [/hr]). An empty bin gets the same ZERO_COUNT_UPPER one-sided 68% Poisson
# upper limit background_subtract() uses instead of sqrt(0) = 0, since this format has one error column and
# no way to say "0 below, 1.14 above". A nan or non-positive exposure gives nan, not a division blow-up
def _final_rate(count, liveTimeSec):
    if not liveTimeSec > 0:
        return float("nan"), float("nan")
    perHour = SECONDS_PER_HOUR / liveTimeSec
    return count * perHour, (np.sqrt(count) if count > 0 else ZERO_COUNT_UPPER) * perHour


# same, background-subtracted with the background scaled to the source exposure. Identical arithmetic to
# background_subtract(), which returns equal low/high errors whenever both counts are non-zero; the
# zero-count branches take the upper limit for the same one-error-column reason as above. NOT floored at
# zero -- the reference file carries negative subtracted rates, and flooring would hide them. A setpoint
# with no background livetime at all comes out nan across every BkgSub column
def _final_bkgsub_rate(count, liveTimeSec, backCount, backTimeSec):
    if not (liveTimeSec > 0 and backTimeSec > 0):
        return float("nan"), float("nan")
    scale = liveTimeSec / backTimeSec
    perHour = SECONDS_PER_HOUR / liveTimeSec
    countErr = np.sqrt(count) if count > 0 else ZERO_COUNT_UPPER
    backErr = (np.sqrt(backCount) if backCount > 0 else ZERO_COUNT_UPPER) * scale
    return (count - backCount * scale) * perHour, np.sqrt(countErr ** 2 + backErr ** 2) * perHour


# per-(p, T), per-rate-column source exposure [s] backed out of the reference file's own numbers, for the
# alternate table. For a plain Poisson rate R = N/tau with error sqrt(N)/tau, R/error^2 = tau, so each of
# the file's five (Rate, Rate Error) pairs carries the exposure that produced it.
#
# Those five do NOT agree inside a row -- the first row implies 1.4171 hr from its Total pair but 1.1614 hr
# from its Wall pair -- so there is no single "the tau" for a setpoint, and taking one per column is the
# only reading that reproduces every pair. The reference file's rates are evidently weighted sums rather
# than count-over-one-exposure, which is what that spread is. Consequence: in the alternate table a row's
# Rate no longer equals the sum of its four region columns, because the five sit on five exposures.
#
# Returns {(p, T): [tau per rate column]}, or None if the file isn't there.
def load_reference_exposures(path):
    if not os.path.exists(path):
        return None
    byPT = {}
    for row in np.atleast_2d(np.genfromtxt(path)):
        taus = []
        for k in range(len(FINAL_RATES_RATE_COLUMNS)):
            rate = row[FINAL_RATES_FIRST_RATE_COLUMN + 2 * k]
            error = row[FINAL_RATES_FIRST_RATE_COLUMN + 2 * k + 1]
            # rate/error^2 is in hours; every exposure below is in seconds
            taus.append(rate / error ** 2 * SECONDS_PER_HOUR
                        if np.isfinite(rate) and np.isfinite(error) and error > 0 else float("nan"))
        byPT[(round(float(row[0]), 4), round(float(row[1]), 4))] = taus
    return byPT


# exposureByPT: {(p, T) rounded to 4dp: [exposure [s] per rate column]} to divide the counts by, or None to
# use each setpoint's own measured livetime for every column. title/extraLines head the summary.txt section
def write_final_rates(savepath, exposureByPT=None, title="Final rate table", extraLines=()):
    if os.path.basename(savepath) == FINAL_RATES_REFERENCE_FILENAME:
        raise ValueError(f"refusing to overwrite the reference file {FINAL_RATES_REFERENCE_FILENAME}")

    sourceByPT = final_rates_source_by_pt()
    columns = ["Pressure [bara]", "Temperature [K]", "Seitz Threshold [keV]"]
    columns += [name for name, _ in FINAL_RATES_MODEL_FIELDS]
    columns += ["Rate [/hr]", "Rate Error [/hr]"]
    columns += [f"{name} Rate{suffix} [/hr]" for name, _ in FINAL_RATES_REGIONS for suffix in ("", " Error")]
    columns += ["BkgSub Rate [/hr]", "BkgSub Rate Error [/hr]"]
    columns += [f"BkgSub {name} Rate{suffix} [/hr]"
                for name, _ in FINAL_RATES_REGIONS for suffix in ("", " Error")]

    rows, missingFields = [], set()
    uncolumnedByRegion = {code: 0 for code in FINAL_RATES_UNCOLUMNED_REGIONS}
    uncolumnedTotal, noExposure = 0, []
    for p, T in pToUse:
        counts, totalCount, liveTimeSec = sourceByPT.get((p, T), ([0] * len(SOURCE_LABELS), 0, 0.0))
        backCounts = backgroundRegionsByPT.get((p, T), [0] * len(SOURCE_LABELS))
        backBins, backTimeSec = backgroundByPT.get((p, T), ([0] * 5, 0.0))
        # the background's true total, the same way: every multiplicity bin, region-independent
        totalBack = sum(backBins)
        for code in uncolumnedByRegion:
            uncolumnedByRegion[code] += counts[code]
        uncolumnedTotal += totalCount - sum(counts[code] for _, code in FINAL_RATES_REGIONS)

        # one exposure per rate column: the measured livetime repeated, or the supplied per-column taus. A
        # setpoint the override has no entry for gets nan across the row rather than silently falling back
        # to the measured livetime, which would mix two exposure conventions in one table
        if exposureByPT is None:
            exposures = [liveTimeSec] * len(FINAL_RATES_RATE_COLUMNS)
        else:
            exposures = exposureByPT.get((round(float(p), 4), round(float(T), 4)))
            if exposures is None:
                noExposure.append((p, T))
                exposures = [float("nan")] * len(FINAL_RATES_RATE_COLUMNS)

        model = seitz_model(p, T)
        row = [float(p), float(T), float(model.Q)]
        for name, attribute in FINAL_RATES_MODEL_FIELDS:
            if attribute is None:
                missingFields.add(name)
            row.append(_final_rates_model_field(model, attribute))

        # Rate is every event at the setpoint, uncut, so it is NOT the sum of the four region columns
        # whenever the combined region is still "other" -- see the module header
        countByColumn = [totalCount] + [counts[code] for _, code in FINAL_RATES_REGIONS]
        backByColumn = [totalBack] + [backCounts[code] for _, code in FINAL_RATES_REGIONS]

        for count, exposure in zip(countByColumn, exposures):
            row += list(_final_rate(count, exposure))
        for count, back, exposure in zip(countByColumn, backByColumn, exposures):
            row += list(_final_bkgsub_rate(count, exposure, back, backTimeSec))
        rows.append(row)

    # default fmt/delimiter/comments reproduce the reference file's encoding exactly: "# " + tab-joined
    # header, then %.18e space-separated
    np.savetxt(savepath, np.asarray(rows, dtype=float), header="\t".join(columns))

    lines = list(extraLines) + [
        f"  wrote {len(rows)} setpoint rows x {len(columns)} columns to {savepath}",
        "  UNCUT: no dome cut, no region exclusion, no multiplicity cut. Rate is every event recorded",
        "  at the setpoint; the region columns break out four of the combined_region() (handscan label",
        "  OR reco position) labels.",
        "  regions: " + ", ".join(f"{name} = combined region {code} ({SOURCE_LABELS[code]})"
                                  for name, code in FINAL_RATES_REGIONS)]
    for code, uncolumned in sorted(uncolumnedByRegion.items()):
        lines.append(f"  combined region {code} ({SOURCE_LABELS[code]}) has no column in this format: "
                     f"{uncolumned} event(s) are in Rate but in none of the region columns")
    unlabelled = uncolumnedTotal - sum(uncolumnedByRegion.values())
    if unlabelled:
        lines.append(f"  {unlabelled} event(s) carry a scan_source code outside {sorted(SOURCE_LABELS)} "
                     f"entirely -- also in Rate, also in no region column")
    lines.append(f"  Rate - (Bulk + Wall + Dome + Bottom) = {uncolumnedTotal} event(s) over all setpoints"
                 + ("" if uncolumnedTotal else "  (so Rate does equal the region sum for this run)"))
    noBackground = [pT for pT in pToUse if backgroundByPT.get(pT, ([0] * 5, 0.0))[1] <= 0]
    if noBackground:
        lines.append(f"  no background livetime at {len(noBackground)} setpoint(s) {noBackground} -- "
                     f"their BkgSub columns are nan")
    if noExposure:
        lines.append(f"  no exposure supplied for {len(noExposure)} setpoint(s) {noExposure} -- "
                     f"their rate columns are nan")
    if missingFields:
        lines.append(f"  {sorted(missingFields)} have no SeitzModel attribute assigned, so those columns "
                     f"are nan -- name the attributes in FINAL_RATES_MODEL_FIELDS to fill them in")
    if rows:
        lines += ["", "  first row, for eyeballing units against the reference file:"]
        lines += [f"    {name:<24}= {value:0.6g}" for name, value in zip(columns, rows[0])]
    report_section(SUMMARY_REPORT, title, lines)


# the primary table: this run's own measured livetime
write_final_rates(final_output_path(FINAL_RATES_FILENAME),
                  title=f"Final rate table ({FINAL_RATES_FILENAME}, measured livetime)")

# and the alternate: same counts, same background, exposures taken from the gray table instead. Skipped
# with a note rather than a crash when that file isn't in GRAY_INPUT_ROOT
_referencePath = os.path.join(GRAY_INPUT_ROOT, FINAL_RATES_REFERENCE_FILENAME)
_referenceExposures = load_reference_exposures(_referencePath)
if _referenceExposures is None:
    report_section(SUMMARY_REPORT, f"Final rate table ({FINAL_RATES_TAU_FILENAME}, tau from the reference file)",
                   [f"  skipped: no {FINAL_RATES_REFERENCE_FILENAME} in {GRAY_INPUT_ROOT}"])
else:
    _tauLines = [f"  exposures backed out of {FINAL_RATES_REFERENCE_FILENAME} as tau = Rate / Rate Error^2,",
                 "  one per rate column, and used in place of this run's measured livetime. The counts,",
                 "  the background and the background livetime are unchanged from the table above.",
                 "",
                 "  tau [hr] per setpoint, per rate column:",
                 "    " + f"{'p':>5} {'T':>6} " + " ".join(f"{n:>8}" for n in FINAL_RATES_RATE_COLUMNS)]
    for _pT in pToUse:
        _taus = _referenceExposures.get((round(float(_pT[0]), 4), round(float(_pT[1]), 4)))
        _shown = ["  n/a  "] * len(FINAL_RATES_RATE_COLUMNS) if _taus is None else \
                 [f"{tau / SECONDS_PER_HOUR:8.4f}" for tau in _taus]
        _tauLines.append(f"    {_pT[0]:5.2f} {_pT[1]:6.1f} " + " ".join(_shown))
    _tauLines += ["",
                  "  those five disagree inside a row, so Rate here does NOT equal the sum of the four",
                  "  region columns -- see load_reference_exposures() for why.",
                  ""]
    write_final_rates(final_output_path(FINAL_RATES_TAU_FILENAME),
                      exposureByPT=_referenceExposures,
                      title=f"Final rate table ({FINAL_RATES_TAU_FILENAME}, tau from the reference file)",
                      extraLines=_tauLines)
## ---------------------------------------------------------------------------------------------------


## ---------------------------------------------------------------------------------------------------
## bubble region fractions: olivia (this script's handscan labels) vs gray (the *_bubble_regions.txt)
##
## The files under GRAY_REGION_DIR carry, per (run date, run index, pressure setpoint), an independent
## reco-based classification of where each bubble was: N Bulk / N Wall / N Dome / N Bottom summing to
## N Localized, with N Above / N 1cam / N Multibub / N Nodet accounting for the rest of N Valid. Both
## identities were checked against the shipped files and hold exactly.
##
## What gets compared, and why it is restricted the way it is:
##
##   * only (run, pressure setpoint) cells present on BOTH sides. Comparing this run's full sample against
##     a gray sample that is missing runs would fold a coverage difference into the fractions.
##   * only single-bubble events on the olivia side, because the gray N Localized sample is single-bubble
##     by construction -- N Multibub is a separate bucket that never enters N Localized.
##   * fractions are taken over Bulk + Wall + Dome + Bottom on both sides. scan_source 4 ("other") has no
##     counterpart inside N Localized, so including it would dilute one denominator and not the other. It
##     is reported separately below, next to the gray Above/1cam/Nodet counts it plausibly overlaps.
##
## So this measures the one thing it can measure cleanly: given a bubble that both analyses localized to
## one of the four regions, do they agree about which region. It is NOT a statement about localization
## efficiency -- that is what the context table at the end of the report is for.
GRAY_REGION_DIR = GRAY_INPUT_ROOT
GRAY_REGION_REPORT = "bubbleRegionFractions.txt"

# gray column name -> the SOURCE_LABELS code the handscan uses for the same place
GRAY_REGION_COLUMNS = [("Bulk", 0), ("Wall", 1), ("Dome", 2), ("Bottom", 3)]
# the gray columns that are not one of the four regions, for the context table
GRAY_OTHER_COLUMNS = ["Above", "1cam", "Multibub", "Nodet"]
# setpoints with fewer localized events than this on either side are kept in the report table but held out
# of the figure. The gray files carry finer setpoints than pToUse (ramp transients at 3.25, 5.75, 6.75 and
# the like), and a cell holding one event plots as a 0% or 100% point with a meaningless error bar that
# rescales the panel and buries every setpoint that actually has statistics
GRAY_REGION_MIN_PLOT_EVENTS = 10
# (this script's side, the txt-file side) -- the names every table heading and plot legend below is built
# from. "olivia" is this pipeline's handscan scan_source labelling; "gray" is the supplied txt files
GRAY_REGION_LABELS = ("olivia", "gray")
# the one-letter column prefixes ("o Bulk" / "g Bulk"), taken from the labels so they cannot drift
GRAY_REGION_PREFIXES = tuple(label[0] for label in GRAY_REGION_LABELS)


# {(run, pset): {column: count}} summed over every gray file, for runs in neutronRuns only. run is the
# "date_index" string the rest of this script uses; pset is rounded to 2dp to key against psetsTemps
def load_gray_region_counts(directory, runList):
    byRunPset = {}
    runSet = set(runList)
    paths = sorted(glob.glob(os.path.join(directory, "*_bubble_regions.txt")))
    filesUsed = set()
    for path in paths:
        with open(path, encoding="utf-8") as handle:
            header = None
            for raw in handle:
                if raw.startswith("#"):
                    header = [c.strip() for c in raw.lstrip("#").strip("\n").split("\t")]
                    continue
                parts = raw.split()
                if len(parts) < len(header):
                    continue
                run = f"{parts[0]}_{parts[1]}"
                if run not in runSet:
                    continue
                filesUsed.add(os.path.basename(path))
                key = (run, round(float(parts[2]), 2))
                cell = byRunPset.setdefault(key, {})
                # header names are "N Bulk", "N Wall", ... -- strip the "N " to key on the bare name
                for name, value in zip(header[3:], parts[3:]):
                    cell[name[2:] if name.startswith("N ") else name] = \
                        cell.get(name[2:] if name.startswith("N ") else name, 0) + int(value)
    return byRunPset, sorted(filesUsed)


# {(run, pset): [count per SOURCE_LABELS code]} for single-bubble handscan events, restricted to the
# (run, pset) cells the gray files also have. Returns the restricted counts and what the restriction cost
def handscan_region_counts(grayKeys):
    byRunPset = {}
    skippedNoGrayCell = 0
    skippedUnstablePset = 0
    for (mult, region), pset, (run, _ev), coord in zip(bubbleCount, psetsTemps, neutronRunEvs,
                                                       sourceSingleBubbleCoords):
        if mult != SINGLE_BUBBLE_MULT:
            continue
        # combined reco + handscan, so this side matches what the wall/dome tables use
        region = combined_region(region, coord)
        if float(pset[0]) != float(pset[1]):
            skippedUnstablePset += 1
            continue
        key = (run, round(float(pset[0]), 2))
        if key not in grayKeys:
            skippedNoGrayCell += 1
            continue
        counts = byRunPset.setdefault(key, [0] * len(SOURCE_LABELS))
        if 0 <= region < len(counts):
            counts[region] += 1
    return byRunPset, skippedUnstablePset, skippedNoGrayCell


# fraction and its binomial error, nan on an empty denominator
def _fraction(count, total):
    if total <= 0:
        return float("nan"), float("nan")
    f = count / total
    return f, np.sqrt(max(f * (1 - f), 0.0) / total)


def _pull(a, aErr, b, bErr):
    denominator = np.sqrt(aErr ** 2 + bErr ** 2)
    if not np.isfinite(denominator) or denominator <= 0:
        return float("nan")
    return (a - b) / denominator


def _pct(value):
    return f"{100 * value:6.2f}" if np.isfinite(value) else "   n/a"


def report_bubble_region_fractions(savepath, figurePath):
    if not os.path.isdir(GRAY_REGION_DIR):
        report_section(SUMMARY_REPORT, "Bubble region fractions vs the gray-region files",
                       [f"  skipped: no {GRAY_REGION_DIR} directory"])
        return

    gray, filesUsed = load_gray_region_counts(GRAY_REGION_DIR, neutronRuns)
    if not gray:
        report_section(SUMMARY_REPORT, "Bubble region fractions vs the gray-region files",
                       [f"  skipped: no gray row matches any run in neutronRuns"])
        return
    hand, skippedUnstable, skippedNoCell = handscan_region_counts(set(gray))

    # a gray cell with no handscan counterpart still gets dropped, so both sides see the same cells
    sharedKeys = sorted(set(gray) & set(hand))
    grayRuns = {run for run, _ in gray}
    missingRuns = sorted(set(neutronRuns) - grayRuns)

    lines = [
        f"{GRAY_REGION_LABELS[0]} = handscan scan_source OR reco position (see combined_region())",
        f"{GRAY_REGION_LABELS[1]} = the reco-based classification in "
        + os.path.basename(GRAY_REGION_DIR),
        f"table columns are prefixed {GRAY_REGION_PREFIXES[0]} = {GRAY_REGION_LABELS[0]} and "
        f"{GRAY_REGION_PREFIXES[1]} = {GRAY_REGION_LABELS[1]}",
        "",
        "restricted to (run, pressure setpoint) cells present in both, to single-bubble events, and to",
        "the four localized regions -- see the block comment in the script for why each restriction is",
        "there. Fractions are over Bulk + Wall + Dome + Bottom on both sides; errors are binomial;",
        f"pulls are ({GRAY_REGION_LABELS[0]} - {GRAY_REGION_LABELS[1]}) / "
        f"sqrt(sigma_{GRAY_REGION_LABELS[0]}^2 + sigma_{GRAY_REGION_LABELS[1]}^2) and treat the two as "
        f"independent,",
        "which they are not (same events, two classifications), so read them as a rough distance.",
        "",
        f"gray files contributing: {len(filesUsed)}",
    ] + [f"    {name}" for name in filesUsed] + [
        "",
        f"  {len(sharedKeys)} (run, setpoint) cells shared, out of {len(gray)} "
        f"{GRAY_REGION_LABELS[1]} and {len(hand)} {GRAY_REGION_LABELS[0]}",
        f"  {len(grayRuns)} of {len(neutronRuns)} neutronRuns runs appear in the gray files",
    ]
    if missingRuns:
        lines.append(f"  {len(missingRuns)} run(s) with no gray row, excluded from both sides: {missingRuns}")
    lines += [f"  {skippedNoCell} single-bubble {GRAY_REGION_LABELS[0]} event(s) fell in no shared cell, "
              f"excluded",
              f"  {skippedUnstable} single-bubble {GRAY_REGION_LABELS[0]} event(s) had "
              f"pset_lo != pset_hi, excluded"]

    # aggregate the shared cells by (pressure setpoint, temperature)
    byPT = {}
    for run, pset in sharedKeys:
        key = (pset, run_temperature(run))
        handTotals, grayTotals, other = byPT.setdefault(key, ([0] * len(SOURCE_LABELS),
                                                              {n: 0 for n, _ in GRAY_REGION_COLUMNS},
                                                              {n: 0 for n in GRAY_OTHER_COLUMNS}))
        for code in range(len(SOURCE_LABELS)):
            handTotals[code] += hand[(run, pset)][code]
        for name, _ in GRAY_REGION_COLUMNS:
            grayTotals[name] += gray[(run, pset)].get(name, 0)
        for name in GRAY_OTHER_COLUMNS:
            other[name] += gray[(run, pset)].get(name, 0)

    regionNames = [name for name, _ in GRAY_REGION_COLUMNS]
    countHeadings = [f"n{label.capitalize()}" for label in GRAY_REGION_LABELS]
    countW = [max(6, len(h)) for h in countHeadings]
    regionBlock = lambda prefix: " ".join(f"{prefix + ' ' + n:>8}" for n in regionNames)
    tableHeader = (f"{'pset':>5} {'T':>6} "
                   + " ".join(f"{h:>{w}}" for h, w in zip(countHeadings, countW)) + " | "
                   + regionBlock(GRAY_REGION_PREFIXES[0]) + " | "
                   + regionBlock(GRAY_REGION_PREFIXES[1]) + " | "
                   + " ".join(f"{'pull ' + n:>9}" for n in regionNames))
    width = len(tableHeader)
    lines += ["", "=" * width,
              "REGION FRACTIONS [%] over Bulk + Wall + Dome + Bottom, by setpoint",
              "=" * width, tableHeader]

    plotRows = []
    for (pset, T) in sorted(byPT):
        handTotals, grayTotals, _other = byPT[(pset, T)]
        handDen = sum(handTotals[code] for _, code in GRAY_REGION_COLUMNS)
        grayDen = sum(grayTotals[name] for name, _ in GRAY_REGION_COLUMNS)
        handF = [_fraction(handTotals[code], handDen) for _, code in GRAY_REGION_COLUMNS]
        grayF = [_fraction(grayTotals[name], grayDen) for name, _ in GRAY_REGION_COLUMNS]
        pulls = [_pull(h, hE, g, gE) for (h, hE), (g, gE) in zip(handF, grayF)]
        lines.append(
            f"{pset:5.2f} {T:6.1f} {handDen:{countW[0]}d} {grayDen:{countW[1]}d} | "
            + " ".join(_pct(f) + "  " for f, _ in handF) + "| "
            + " ".join(_pct(f) + "  " for f, _ in grayF) + "| "
            + " ".join(f"{p:9.2f}" if np.isfinite(p) else "      n/a" for p in pulls))
        plotRows.append({"pset": pset, "T": T, "hand": handF, "gray": grayF,
                         "handDen": handDen, "grayDen": grayDen})

    # pooled, per temperature and overall -- the headline numbers
    lines += ["", "=" * width, "POOLED over every shared cell", "=" * width,
              f"{'sample':>12} "
              + " ".join(f"{h:>{w}}" for h, w in zip(countHeadings, countW)) + " | "
              + regionBlock(GRAY_REGION_PREFIXES[0]) + " | "
              + regionBlock(GRAY_REGION_PREFIXES[1]) + " | "
              + " ".join(f"{'pull ' + n:>9}" for n in regionNames)]
    temperatures = sorted({T for _, T in byPT})
    for label, keys in [(f"{T:0.1f} K", [k for k in byPT if k[1] == T]) for T in temperatures] \
            + [("all", list(byPT))]:
        handTotals = [sum(byPT[k][0][code] for k in keys) for code in range(len(SOURCE_LABELS))]
        grayTotals = {name: sum(byPT[k][1][name] for k in keys) for name, _ in GRAY_REGION_COLUMNS}
        handDen = sum(handTotals[code] for _, code in GRAY_REGION_COLUMNS)
        grayDen = sum(grayTotals.values())
        handF = [_fraction(handTotals[code], handDen) for _, code in GRAY_REGION_COLUMNS]
        grayF = [_fraction(grayTotals[name], grayDen) for name, _ in GRAY_REGION_COLUMNS]
        pulls = [_pull(h, hE, g, gE) for (h, hE), (g, gE) in zip(handF, grayF)]
        lines.append(
            f"{label:>12} {handDen:{countW[0]}d} {grayDen:{countW[1]}d} | "
            + " ".join(_pct(f) + "  " for f, _ in handF) + "| "
            + " ".join(_pct(f) + "  " for f, _ in grayF) + "| "
            + " ".join(f"{p:9.2f}" if np.isfinite(p) else "      n/a" for p in pulls))

    # context: what each side does with the events that are NOT one of the four localized regions. These
    # denominators are different things, so this is orientation, not a comparison
    handOther = sum(sum(byPT[k][0][code] for code in range(len(SOURCE_LABELS))
                        if code not in {c for _, c in GRAY_REGION_COLUMNS}) for k in byPT)
    handAll = sum(sum(byPT[k][0]) for k in byPT)
    grayOther = {name: sum(byPT[k][2][name] for k in byPT) for name in GRAY_OTHER_COLUMNS}
    grayLocalized = sum(sum(byPT[k][1].values()) for k in byPT)
    grayValid = grayLocalized + sum(grayOther.values())
    lines += ["", "=" * width,
              "CONTEXT -- events outside the four localized regions (different denominators, not a comparison)",
              "=" * width,
              f"  {GRAY_REGION_LABELS[0]}: {handOther} of {handAll} single-bubble events are scan_source "
              f"{sorted(set(SOURCE_LABELS) - {c for _, c in GRAY_REGION_COLUMNS})} "
              f"= {100 * handOther / handAll if handAll else float('nan'):0.2f}% "
              f"(these are excluded from the fractions above)",
              f"  {GRAY_REGION_LABELS[1]}: {grayLocalized} of {grayValid} valid events are localized "
              f"= {100 * grayLocalized / grayValid if grayValid else float('nan'):0.2f}%"]
    for name in GRAY_OTHER_COLUMNS:
        lines.append(f"      N {name:<9}= {grayOther[name]:5d}  "
                     f"({100 * grayOther[name] / grayValid if grayValid else float('nan'):5.2f}% of valid)")
    lines.append(f"  note: the {GRAY_REGION_LABELS[0]} side above is single-bubble only, so "
                 f"{GRAY_REGION_LABELS[1]} N Multibub has no {GRAY_REGION_LABELS[0]} counterpart here "
                 f"by construction")

    thin = [(r["pset"], r["T"], min(r["handDen"], r["grayDen"])) for r in plotRows
            if min(r["handDen"], r["grayDen"]) < GRAY_REGION_MIN_PLOT_EVENTS]
    if thin:
        lines += ["", f"  {len(thin)} setpoint(s) have < {GRAY_REGION_MIN_PLOT_EVENTS} localized events on "
                      f"one side and are in the table above but held out of the figure:",
                  "    " + ", ".join(f"{p:0.2f} bara / {T:0.1f} K (n={n})" for p, T, n in thin)]

    with open(savepath, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    report_section(SUMMARY_REPORT, "Bubble region fractions vs the gray-region files",
                   [f"  wrote {savepath}",
                    f"  {len(sharedKeys)} shared (run, setpoint) cells, "
                    f"{sum(r['handDen'] for r in plotRows)} {GRAY_REGION_LABELS[0]} / "
                    f"{sum(r['grayDen'] for r in plotRows)} {GRAY_REGION_LABELS[1]} "
                    f"single-bubble localized events"]
                   # the pooled block: its column header plus one row per temperature plus "all"
                   + lines[lines.index("POOLED over every shared cell") + 2:][:len(temperatures) + 2])
    plot_bubble_region_fractions(plotRows, regionNames, figurePath)


# one panel per (temperature, region): fraction vs pressure setpoint, olivia against gray
def plot_bubble_region_fractions(plotRows, regionNames, savepath):
    plotRows = [r for r in plotRows
                if min(r["handDen"], r["grayDen"]) >= GRAY_REGION_MIN_PLOT_EVENTS]
    if not plotRows:
        return
    temperatures = sorted({r["T"] for r in plotRows})
    fig, axes = plt.subplots(len(temperatures), len(regionNames),
                             figsize=(4.0 * len(regionNames), 3.2 * len(temperatures)), squeeze=False)
    for rowIdx, T in enumerate(temperatures):
        rows = sorted([r for r in plotRows if r["T"] == T], key=lambda r: r["pset"])
        psets = [r["pset"] for r in rows]
        for colIdx, name in enumerate(regionNames):
            ax = axes[rowIdx][colIdx]
            ax.errorbar(psets, [100 * r["hand"][colIdx][0] for r in rows],
                        yerr=[100 * r["hand"][colIdx][1] for r in rows],
                        fmt="o", color="black", markersize=4, linewidth=1, capsize=2,
                        label=GRAY_REGION_LABELS[0])
            ax.errorbar([p + 0.03 for p in psets], [100 * r["gray"][colIdx][0] for r in rows],
                        yerr=[100 * r["gray"][colIdx][1] for r in rows],
                        fmt="s", color="tab:red", markersize=4, linewidth=1, capsize=2,
                        label=GRAY_REGION_LABELS[1])
            ax.set_title(f"{name}, {T:0.1f} K", fontsize=9)
            ax.set_xlabel("pressure setpoint [bara]", fontsize=8)
            ax.set_ylabel("fraction of localized [%]", fontsize=8)
            ax.tick_params(labelsize=7)
            ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(savepath, dpi=130)
    plt.close(fig)


report_bubble_region_fractions(output_path("comparison", GRAY_REGION_REPORT),
                               output_path("comparison", "bubbleRegionFractions.png"))
## ---------------------------------------------------------------------------------------------------


## ---------------------------------------------------------------------------------------------------
## how much of the sample is wall or dome
##
## Same counts the rate tables above are built from -- uncut, single- and multi-bubble alike -- so these
## percentages describe those tables directly rather than a separately selected sample. Wall and dome are
## the two surface-adjacent regions, so their combined share is the number worth watching: it is the part
## of the rate the dome cut and any future wall cut would be acting on.
##
## The denominator is the Rate column's own total: every event at the setpoint, including the "other"
## scan_source events that have no region column. That makes "% wall or dome" read directly against the
## Rate/Wall Rate/Dome Rate columns of the tables in output/. It is NOT the fraction of localized events --
## report_bubble_region_fractions() below does that, over Bulk+Wall+Dome+Bottom, for the gray comparison.
WALL_OR_DOME_FILENAME = "wallOrDomeFractions.txt"
WALL_OR_DOME_REGIONS = [("Wall", 1), ("Dome", 2)]


# count -> (percent, binomial error on the percent), nan on an empty denominator
def _percent(count, total):
    if total <= 0:
        return float("nan"), float("nan")
    f = count / total
    return 100 * f, 100 * np.sqrt(max(f * (1 - f), 0.0) / total)


# width of a "value +- error" percentage cell, so the headings and the numbers under them agree
PCT_CELL_WIDTH = 15


def _pct_pm(count, total):
    value, error = _percent(count, total)
    cell = f"{value:6.2f} +-{error:5.2f}" if np.isfinite(value) else "n/a"
    return cell.rjust(PCT_CELL_WIDTH)


# the four samples the cross-comparison table puts side by side, each as {temperature: per-region counts}.
# olivia's two come from the handscan scan_source labels this pipeline already loads; gray's two come from
# the *_bubble_regions.txt files, which cover the background runs as well as the source runs
def _wall_or_dome_samples():
    samples = []
    setpoints = set(pToUse)

    # SINGLE-BUBBLE ONLY on this side, to match gray's N Localized. The all-multiplicity counts the rate
    # tables use are deliberately not reused here: the handscan gives one scan_source per event, so a
    # multi-bubble event carries a single region label, and gray never puts such an event in a region
    # column at all -- it goes in N Multibub. Counting them on one side only would bias every share
    byT = {}
    for (p, T), counts in sourceCombinedRegionsByPT.items():
        if (p, T) not in setpoints:
            continue
        acc = byT.setdefault(T, [0] * len(SOURCE_LABELS))
        for code, count in enumerate(counts):
            acc[code] += count
    samples.append((f"{GRAY_REGION_LABELS[0]} source", byT))

    byT = {}
    for (_p, T), counts in backgroundCombinedRegionsByPT.items():
        acc = byT.setdefault(T, [0] * len(SOURCE_LABELS))
        for code, count in enumerate(counts):
            acc[code] += count
    samples.append((f"{GRAY_REGION_LABELS[0]} background", byT))

    for label, runList in ((f"{GRAY_REGION_LABELS[1]} source", neutronRuns),
                           (f"{GRAY_REGION_LABELS[1]} background", backgroundList)):
        cells, _filesUsed = load_gray_region_counts(GRAY_REGION_DIR, runList)
        byT = {}
        for (run, _pset), cell in cells.items():
            acc = byT.setdefault(run_temperature(run), [0] * len(SOURCE_LABELS))
            for name, code in FINAL_RATES_REGIONS:
                acc[code] += cell.get(name, 0)
        samples.append((label, byT))
    return samples


# wall / dome / wall-or-dome as a share of the four localized regions, for source and background on both
# sides. The denominator is Bulk + Wall + Dome + Bottom rather than the Rate total used above, because
# that is the one denominator both sides can form identically: gray's N Localized IS that sum, and gray
# has no counterpart for the handscan's "other" bucket. So the numbers here are NOT the same percentages
# as the first table -- they answer "of the bubbles we placed, how many were surface", for four samples.
# Single-bubble only on both sides, for the same reason: gray's N Localized excludes multi-bubble events.
#
# The background rows are the point of the table. A wall or dome share that is the same source-on as
# source-off says those events are not neutron-induced; a source-on excess says they are.
def _wall_or_dome_cross_comparison(samples):
    names = [name for name, _ in WALL_OR_DOME_REGIONS]
    header = (f"{'sample':>20} {'T':>7} {'N loc':>7} "
              + " ".join(f"{'N ' + n:>6}" for n in names) + f" {'N either':>8} | "
              + " ".join(f"{'% ' + n:>{PCT_CELL_WIDTH}}" for n in names)
              + f" | {'% wall or dome':>{PCT_CELL_WIDTH}}")
    lines = ["=" * len(header),
             "WALL / DOME SHARE OF LOCALIZED EVENTS -- source vs background, olivia vs gray",
             "=" * len(header),
             "single-bubble only, over Bulk + Wall + Dome + Bottom. Different denominator AND different",
             "sample from the table above, which is all multiplicities over the Rate total -- these two",
             f"tables are not meant to agree. {GRAY_REGION_LABELS[0]} wall/dome are the handscanner's",
             f"scan_source labels; {GRAY_REGION_LABELS[1]}'s are reco-based. No geometric dome cut here.",
             "=" * len(header), header, "-" * len(header)]
    for label, byT in samples:
        allCounts = [0] * len(SOURCE_LABELS)
        for T in sorted(byT) + ["all"]:
            counts = allCounts if T == "all" else byT[T]
            if T != "all":
                allCounts = [a + b for a, b in zip(allCounts, counts)]
            localized = sum(counts[code] for _name, code in FINAL_RATES_REGIONS)
            regionCounts = [counts[code] for _name, code in WALL_OR_DOME_REGIONS]
            both = sum(regionCounts)
            lines.append(
                f"{label if T == sorted(byT)[0] else '':>20} "
                f"{'all' if T == 'all' else f'{T:0.1f} K':>7} {localized:7d} "
                + " ".join(f"{c:6d}" for c in regionCounts) + f" {both:8d} | "
                + " ".join(_pct_pm(c, localized) for c in regionCounts) + " | "
                + _pct_pm(both, localized))
        lines.append("-" * len(header))
    return lines


def report_wall_or_dome_fractions(savepath):
    sourceByPT = {pT: (counts, sum(counts), 0.0) for pT, counts in sourceCombinedRegionsByPT.items()}
    names = [name for name, _ in WALL_OR_DOME_REGIONS]
    header = (f"{'pset':>5} {'T':>6} {'N':>6} "
              + " ".join(f"{'N ' + n:>6}" for n in names) + f" {'N either':>8} | "
              + " ".join(f"{'% ' + n:>{PCT_CELL_WIDTH}}" for n in names)
              + f" | {'% wall or dome':>{PCT_CELL_WIDTH}}")
    diag = sourceCombinedDiagnostics
    lines = [
        "percentage of single-bubble events that are wall or dome",
        "",
        "region = handscan scan_source OR reco position, whichever calls it wall/dome first:",
        f"  dome   scan_source == {DOME_REGION_CODE}, or reco z > {DOME_Z_THRESHOLD_CM * 10:0.1f} mm",
        f"  wall   scan_source == {WALL_REGION_CODE}, or reco radius >= {WALL_RADIUS_MM:0.3f} mm",
        f"         ({WALL_PROXIMITY_MM:0.0f} mm inside the {WALL_INNER_RADIUS_MM:0.3f} mm inner surface "
        f"and outward, so anything",
        f"         reconstructed inside the {WALL_OUTER_RADIUS_MM - WALL_INNER_RADIUS_MM:0.2f} mm of glass "
        f"is included), for {CYLINDER_BOTTOM_MM:0.1f} <= z <= {CYLINDER_TOP_MM:0.1f} mm",
        "  dome wins where both fire -- the dome threshold reaches 13.8 mm below the cylinder top, so the",
        "  top of the wall is inside it. The overlap is counted below rather than hidden.",
        "",
        "SINGLE-BUBBLE ONLY. The handscan labels a multi-bubble event as a whole and such events have no",
        "reliable 3D coord, so neither input means anything for them. That makes N here smaller than the",
        "Rate column of the tables in " + os.path.basename(OUTPUT_ROOT) + "/, which counts every event.",
        "",
        f"  {diag['single']} single-bubble events at a fixed setpoint",
        f"  {diag['noCoord']} of them have no valid reco coord -- handscan label only for those",
        f"  {diag['bothFlagged']} satisfy the wall AND dome tests at once; dome wins for those, so they",
        "  appear in the dome column, not the wall one (this is the OVERLAP -- the 'N either' column",
        "  below is the UNION, and the two regions are exclusive there)",
        f"  {diag['movedByReco']} were moved out of their handscan label by the reco position",
        "",
        "Errors are binomial. No dome cut is applied: this describes what the uncut sample contains,",
        "which is what makes it useful for sizing a cut.",
        "",
        "=" * len(header), header, "=" * len(header),
    ]

    pooled = {}
    for pT in pToUse:
        counts, totalCount, _liveTimeSec = sourceByPT.get(pT, ([0] * len(SOURCE_LABELS), 0, 0.0))
        regionCounts = [counts[code] for _, code in WALL_OR_DOME_REGIONS]
        both = sum(regionCounts)
        lines.append(
            f"{pT[0]:5.2f} {pT[1]:6.1f} {totalCount:6d} "
            + " ".join(f"{c:6d}" for c in regionCounts) + f" {both:8d} | "
            + " ".join(_pct_pm(c, totalCount) for c in regionCounts) + " | "
            + _pct_pm(both, totalCount))
        key = pT[1]
        acc = pooled.setdefault(key, [0] * (len(names) + 1))
        for i, c in enumerate(regionCounts):
            acc[i] += c
        acc[-1] += totalCount

    lines += ["", "=" * len(header), "POOLED", "=" * len(header),
              f"{'sample':>12} {'N':>6} " + " ".join(f"{'N ' + n:>6}" for n in names)
              + f" {'N either':>8} | " + " ".join(f"{'% ' + n:>{PCT_CELL_WIDTH}}" for n in names)
              + f" | {'% wall or dome':>{PCT_CELL_WIDTH}}"]
    allAcc = [0] * (len(names) + 1)
    for T in sorted(pooled) + ["all"]:
        acc = allAcc if T == "all" else pooled[T]
        if T != "all":
            allAcc = [a + b for a, b in zip(allAcc, acc)]
        total, regionCounts = acc[-1], acc[:-1]
        both = sum(regionCounts)
        label = "all" if T == "all" else f"{T:0.1f} K"
        lines.append(
            f"{label:>12} {total:6d} " + " ".join(f"{c:6d}" for c in regionCounts) + f" {both:8d} | "
            + " ".join(_pct_pm(c, total) for c in regionCounts) + " | "
            + _pct_pm(both, total))

    lines += ["", ""] + _wall_or_dome_cross_comparison(_wall_or_dome_samples())

    with open(savepath, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    # the pooled block is the headline, so it goes in summary.txt too rather than only in its own file
    report_section(SUMMARY_REPORT, "Percentage of events that are wall or dome",
                   [f"  wrote {savepath}", ""] + lines[-(len(pooled) + 4):])


report_wall_or_dome_fractions(final_output_path(WALL_OR_DOME_FILENAME))
## ---------------------------------------------------------------------------------------------------



# how the (p, T) pairs get partitioned into plotted groups, as (label, [(p, T), ...]) pairs. Normally one
# group per pair; with useAlternateSeitzBinning on, every pair landing in the same alternateSeitzBinEdges
# range is pooled into one group, and pairs outside the edges (in practice the high-Seitz ones above the
# top edge) each keep their own group. The label goes into every per-group filename
def build_pt_bins():
    if not useAlternateSeitzBinning:
        return [(f"{p}{T}", [(p, T)]) for p, T in pToUse]

    lowest, highest = alternateSeitzBinEdges[0], alternateSeitzBinEdges[-1]
    bins = []
    for lo, hi in zip(alternateSeitzBinEdges, alternateSeitzBinEdges[1:]):
        members = [pT for pT in pToUse if lo <= seitz_threshold(*pT) < hi]
        if members:
            bins.append((f"{lo:g}-{hi:g}keV", members))
    bins += [(f"{p}{T}", [(p, T)]) for p, T in pToUse
             if not (lowest <= seitz_threshold(p, T) < highest)]
    return bins


# Seitz threshold and sim-predicted multiplicity counts for one bin of (p, T) pairs. A lone pair is just its
# own Seitz Q with the sim evaluated there -- identical to what build_groups() did before. A pooled bin
# averages both the threshold and the per-member sim counts weighted by each member's source livetime,
# since the pooled data rate is itself the livetime-weighted average of its members' rates
def pooled_seitz_prediction(members, simExcludedRegions):
    seitzByMember = [seitz_threshold(p, T) for p, T in members]
    countsByMember = [get_multiplicity_counts(seitz, excludedRegionsOverride=simExcludedRegions)[0]
                      for seitz in seitzByMember]

    # EQUAL weights -- every member setpoint counts the same regardless of how long it ran. The
    # livetime-weighted version this replaces is in neutronExcluded.py; see the header for what changes
    n = len(members)
    seitz = sum(seitzByMember) / n
    # get_multiplicity_counts() can come back shorter than multiplicity_cut when a member's sim tops out
    # at a low multiplicity, so pad the short ones with zeros before averaging
    numBins = max(len(counts) for counts in countsByMember)
    seitzCounts = [
        sum(counts[i] if i < len(counts) else 0.0 for counts in countsByMember) / n
        for i in range(numBins)
    ]
    return seitz, seitzCounts


# "P = 1.50 bar, T = 116.7 K" for a normal one-(p,T) group, or the pooled members for an alternate-binning group
def group_pt_description(g):
    if len(g["members"]) == 1:
        p, T = g["members"][0]
        return f"P = {p:0.2f} bar, T = {T:0.1f} K"
    return f"{g['label']} bin: " + ", ".join(f"({p:0.2f} bar, {T:0.1f} K)" for p, T in g["members"])


# short "which group is this" tag for the chi2 debug dump (the avg-group pipeline passes a group with no members)
def group_debug_label(dataGroup):
    if "members" not in dataGroup:
        return "(avg group)"
    if len(dataGroup["members"]) == 1:
        return f'(p={dataGroup["p"]:0.2f}, T={dataGroup["T"]:0.1f})'
    return f'(pooled {dataGroup["label"]}: {len(dataGroup["members"])} (p, T) pairs)'


# subtract + floor + rate-convert one set of multiplicity bins. Called once on the 5 raw classes and once on
# the 3 rebinned ones -- NOT a quadrature sum of the 5-bin errors, because the zero-count upper limit is a
# property of an empty bin, and quadrature-summing three of them would hand a 3+ bin holding one real count
# an upper error of sqrt(1 + 1.14^2 + 1.14^2) instead of the sqrt(1) its own count actually implies
def _rate_fields(binCounts, backgroundBinCounts, liveTimeSec, backgroundTime, suffix):
    (backBins, backErrorLow, backErrorHigh, backSubBins, backSubErrorLow, backSubErrorHigh,
     binCountErrorLow, binCountErrorHigh) = \
        background_subtract(binCounts, liveTimeSec, backgroundBinCounts, backgroundTime)

    flooredBins, flooredErrorLow, flooredErrorHigh = floor_at_zero(backSubBins, backSubErrorLow, backSubErrorHigh)

    liveTimeMin = liveTimeSec / 60
    toRate = lambda values: [v / liveTimeMin for v in values]
    return {
        "binCounts" + suffix: toRate(binCounts),
        "binCountErrorLow" + suffix: toRate(binCountErrorLow),
        "binCountErrorHigh" + suffix: toRate(binCountErrorHigh),
        "backBins" + suffix: toRate(backBins),
        "backErrorLow" + suffix: toRate(backErrorLow),
        "backErrorHigh" + suffix: toRate(backErrorHigh),
        # plotted (floored at 0) ...
        "backSub" + suffix: toRate(flooredBins),
        "backSubErrorLow" + suffix: toRate(flooredErrorLow),
        "backSubErrorHigh" + suffix: toRate(flooredErrorHigh),
        # ... and the same rates as actually subtracted, which the fits and chi2 use
        "backSubUnfloored" + suffix: toRate(backSubBins),
        "errLowUnfloored" + suffix: toRate(backSubErrorLow),
        "errHighUnfloored" + suffix: toRate(backSubErrorHigh),
    }


# bins the exposures inGroup() selects, vetoing the bubbles isVetoed() rejects but keeping their live time,
# background-subtracts, and converts to a rate in counts/minute -- shared core of both build_groups()'s
# per-(p,T) loop and the avg-group section. backgroundIsVetoed is isVetoed's counterpart for the
# background sample -- the same dome exclusion, applied the same way (bubble gone, live time kept), so the
# background being subtracted describes the same region the source count does
def compute_rate_group(inGroup, isVetoed, backgroundPTs, backgroundIsVetoed):
    binCounts, liveTimeSec = bin_multiplicities(bubbleCount, sourceTimes, inGroup, isVetoed)
    backgroundBinCounts, backgroundTime = background_for(backgroundPTs, backgroundIsVetoed)
    return {
        "liveTime": liveTimeSec / 60,
        "liveTimeSec": liveTimeSec,
        # raw (un-scaled, un-rated) inputs to the subtraction, for the chi2 debug dump
        "binCountsRaw": binCounts,
        "backCountsRaw": list(backgroundBinCounts),
        "backgroundTimeSec": backgroundTime,
        "backgroundScale": liveTimeSec / backgroundTime if backgroundTime > 0 else 0.0,
        # 5 multiplicity classes (1,2,3,4,5+) ...
        **_rate_fields(binCounts, backgroundBinCounts, liveTimeSec, backgroundTime, "Full"),
        # ... and the 3 plotted ones (1,2,3+), re-derived from the rebinned raw counts
        **_rate_fields(rebin(binCounts), rebin(backgroundBinCounts), liveTimeSec, backgroundTime, ""),
    }


# builds the full per-(p,T) "groups" list for one (simExcludedRegions, positionDomeFlags) combination -- lets
# the pre/post pipeline independently control the sim-side cut and real-side position cut.
# backgroundPositionDomeFlags is positionDomeFlags's counterpart for the background sample (built by
# compute_position_dome_flags() against backgroundBubbleCount/backgroundRunEvs rather than the source's
# own arrays) -- the caller passes the matching one, so the background is dome-cut the same way the source
# is for this pipeline, not left uncut while the source isn't
def build_groups(simExcludedRegions, positionDomeFlags, backgroundPositionDomeFlags):
    groups = []
    is_region_excluded = make_is_region_excluded(positionDomeFlags)
    background_is_region_excluded = make_is_region_excluded(backgroundPositionDomeFlags)

    for label, members in build_pt_bins():
        memberSet = set(members)
        rateGroup = compute_rate_group(
            inGroup=lambda i, memberSet=memberSet: (psetsTemps[i][0], psetsTemps[i][2]) in memberSet,
            isVetoed=is_region_excluded,
            backgroundPTs=members,
            backgroundIsVetoed=background_is_region_excluded,
        )

        # seitz threshold for this group, fed straight into the dome-excluded Cf sim counts
        seitz, seitzCounts = pooled_seitz_prediction(members, simExcludedRegions)
        groups.append({
            "label": label,
            "members": members,
            # p/T only exist for a plain one-setpoint group -- a pooled bin has several of each
            **({"p": members[0][0], "T": members[0][1]} if len(members) == 1 else {}),
            "seitz": seitz,
            "seitzCountsFull": seitzCounts,
            "errLow": rateGroup["backSubErrorLow"],
            "errHigh": rateGroup["backSubErrorHigh"],
            "seitzRate": rebin(seitzCounts),
            "seitzCounts": rebin(seitzCounts),
            "bestFit": {},
            **rateGroup,
        })

    # sort from lowest to highest seitz threshold
    groups.sort(key=lambda g: g["seitz"])

    # FIT normalizationFactor to the measured rates, livetime-weighted per group so long, low-noise
    # exposures pull the shared factor harder than short ones -- see the module header for the formula
    # and why this is the alternative to activityNormalizationFactor. g["backSub"] and g["seitzCounts"]
    # are both already rebinned (3 bins: 1, 2, 3+), so ratioByGroup is over the same bins in both
    observedRatesByGroup = [g["backSub"] for g in groups]
    simCountsByGroup = [g["seitzCounts"] for g in groups]
    liveTimeByGroup = [g["liveTime"] for g in groups]
    ratioByGroup = [
        sum(observed) / sum(sim) for observed, sim in zip(observedRatesByGroup, simCountsByGroup)
    ]
    normalizationFactor = (
        sum(t * ratio for t, ratio in zip(liveTimeByGroup, ratioByGroup)) / sum(liveTimeByGroup)
    )
    for g in groups:
        g["seitzRate"] = [normalizationFactor * ratio for ratio in g["seitzCounts"]]

    return groups, normalizationFactor


# "withoutDomeCut": sim cut off, real side at the loose threshold. "withDomeCut": sim cut on, real side at the tighter threshold. Each gets its own run_pipeline() subfolder
groupsWithoutDomeCut, normalizationFactorWithoutDomeCut = build_groups(
    [], positionDomeFlagsPre, backgroundPositionDomeFlagsPre)
groupsWithDomeCut, normalizationFactorWithDomeCut = build_groups(
    ["dome"], positionDomeFlagsPost, backgroundPositionDomeFlagsPost)

# combinedMultiplicityFinal.pdf drops the two highest-threshold groups (1.91 and 3.10 keV); groups are
# sorted by seitz, so that is the tail. Kept in one place so the figure and the report describing it can
# never disagree about which panels are in play
FINAL_PLOT_SLICE = slice(0, -2)

# fully uncut baseline: sim off, real side has neither the PRE (always-applied) nor the POST dome
# threshold -- combined-plot-only, doesn't get the full run_pipeline() treatment (no linhists/scans/Z-dists)
positionDomeFlagsNone = [False] * len(bubbleCount)
backgroundPositionDomeFlagsNone = [False] * len(backgroundBubbleCount)
groupsNoCutAtAll, normalizationFactorNoCutAtAll = build_groups(
    [], positionDomeFlagsNone, backgroundPositionDomeFlagsNone)

_decayYears = (ACTIVITY_EVALUATION_DATE - CF252_REFERENCE_DATE).days / 365.25
report_section(SUMMARY_REPORT, "Source activity, for reference (this variant does NOT use it directly)", [
    f"  {CF252_REFERENCE_DATE} -> {ACTIVITY_EVALUATION_DATE} = "
    f"{(ACTIVITY_EVALUATION_DATE - CF252_REFERENCE_DATE).days:,} days = {_decayYears:0.4f} yr",
    f"  {_decayYears:0.4f} yr / {CF252_HALF_LIFE_YEARS} yr = "
    f"{_decayYears / CF252_HALF_LIFE_YEARS:0.4f} half-lives  ->  "
    f"decay factor {2 ** (-_decayYears / CF252_HALF_LIFE_YEARS):0.6f}",
    f"  {CF252_NEUTRONS_PER_SECOND_AT_REFERENCE:,.0f} n/s at reference  ->  "
    f"{cf252_neutrons_per_minute() / 60:0.2f} n/s  ->  {cf252_neutrons_per_minute():,.1f} n/min today",
    f"  {SIMULATED_NEUTRONS:,.0f} simulated neutrons / {cf252_neutrons_per_minute():,.1f} n/min = "
    f"{simulatedMinutes:,.1f} min of equivalent running",
    f"  activity-normalized factor = 1/{simulatedMinutes:,.1f} = {activityNormalizationFactor:0.6g} "
    f"count/min per sim count",
    "",
    "  neutronExcluded.py applies that factor directly, with no fit. THIS variant instead fits",
    "  normalizationFactor to the measured rates -- livetime-weighted per pipeline, so long exposures",
    "  pull the shared factor harder than short ones -- and only prints the activity number above so it",
    "  can be compared against what the fit below actually lands on. Run both scripts and compare",
    "  combinedMultiplicityFinal.pdf to see what the fit is absorbing that the activity number doesn't.",
    "",
    "  Pooling inside a Seitz bin also uses EQUAL member weights here, not livetime. The data side of a",
    "  pooled bin is still exposure-weighted by construction, so the two sides are combined differently --",
    "  see the module header.",
    "",
    "  activity number assumes the 50M neutrons are generated into 4pi at the source, and that Table I's",
    "  source is this source. Sensitivity to ACTIVITY_EVALUATION_DATE is ~2%/month -- see the constant's",
    "  comment. None of that matters to the fitted factor below, which never reads this number."])

report_section(SUMMARY_REPORT, "Sim-to-data normalization factor weighting the simulated data",
               [f"    fitted, livetime-weighted per pipeline (see the module header for the formula);",
                f"    units = count/min of data per raw sim count. For reference, the activity-based",
                f"    factor above is {activityNormalizationFactor:0.6g}"]
               + [f"  {label:<34}: {factor:0.6g}" for label, factor in
                  [("withoutDomeCut (sim cut off)", normalizationFactorWithoutDomeCut),
                   ("withDomeCut (sim cut on)", normalizationFactorWithDomeCut),
                   ("noCutAtAll (combined plots only)", normalizationFactorNoCutAtAll)]])

# every group's background-subtracted data rate and Seitz/sim-predicted rate, pre-cut vs post-cut
def report_rates(groupsPre, groupsPost):
    binLabels = rebinnedBinLabels
    lines = ["data = background-subtracted real rate, sim = Seitz/sim-predicted rate",
             "the two livetimes are equal by construction: the cut removes dome BUBBLES, not the "
             "exposures they came from, so the denominator is the same on both sides", ""]
    for gPre, gPost in zip(groupsPre, groupsPost):
        lines.append(f"Seitz = {gPre['seitz']:0.2f} keV  ({group_pt_description(gPre)})")
        lines.append(f"  Pre-cut  (Z > {DOME_Z_THRESHOLD_CM_PRE}cm, sim off):  livetime = {gPre['liveTime']:0.2f} min")
        for label, rate, errLow, errHigh, simRate in zip(
                binLabels, gPre["backSub"], gPre["errLow"], gPre["errHigh"], gPre["seitzRate"]):
            lines.append(f"    mult {label:<3}: data = {rate:0.4f} (+{errHigh:0.4f}/-{errLow:0.4f})   sim = {simRate:0.4f}")
        lines.append(f"  Post-cut (Z > {DOME_Z_THRESHOLD_CM_POST}cm, sim on): livetime = {gPost['liveTime']:0.2f} min")
        for label, rate, errLow, errHigh, simRate in zip(
                binLabels, gPost["backSub"], gPost["errLow"], gPost["errHigh"], gPost["seitzRate"]):
            lines.append(f"    mult {label:<3}: data = {rate:0.4f} (+{errHigh:0.4f}/-{errLow:0.4f})   sim = {simRate:0.4f}")
        lines.append("")
    report_section(GROUP_REPORT,
                   f"Rates pre-cut (uncut baseline: real side Z > {DOME_Z_THRESHOLD_CM_PRE}cm, sim off) "
                   f"vs post-cut (real side Z > {DOME_Z_THRESHOLD_CM_POST}cm, sim on), count/min",
                   lines)

report_rates(groupsWithoutDomeCut, groupsWithDomeCut)


# diagnostic for the zero-count branch of background_subtract(): a bin with no source-on events is drawn at
# 0 rate with only the ZERO_COUNT_UPPER limit above it, so this lists every such bin next to the raw
# source-on / source-off event counts and livetimes behind it -- enough to confirm the 0 is a genuinely
# empty bin and not a livetime, (p, T) matching or rebinning slip
def report_empty_bins(groupLists):
    fullLabels = ["1", "2", "3", "4", "5+"]
    lines = ["a bin with no source-on events is drawn at 0 rate with only the "
             f"{ZERO_COUNT_UPPER:0.4f}-count Poisson upper limit above it"]
    for pipelineLabel, groups in groupLists:
        lines += ["", f"=== {pipelineLabel} ==="]
        emptyBins, allBins = 0, 0
        for g in groups:
            sourceCounts, backCounts = g["binCountsRaw"], g["backCountsRaw"]
            scale = g["backgroundScale"]
            # the 5 raw multiplicity classes, plus the rebinned 3+ the plots and ratesPrePostCut.txt draw
            bins = list(zip(fullLabels, sourceCounts, backCounts))
            bins.append((rebinnedBinLabels[-1], rebin(sourceCounts)[-1], rebin(backCounts)[-1]))
            allBins += len(bins)
            empty = [(label, n, b) for label, n, b in bins if n == 0]
            if not empty:
                continue
            emptyBins += len(empty)
            sourceEvents, backEvents = sum(sourceCounts), sum(backCounts)
            lines.append(f"  Seitz = {g['seitz']:0.2f} keV  ({group_pt_description(g)})")
            lines.append(f"    livetime: source-on = {g['liveTimeSec']:0.1f} s, "
                         f"source-off = {g['backgroundTimeSec']:0.1f} s  (background scale = {scale:0.4f})")
            lines.append(f"    group totals over mult 1..5+: source-on = {sourceEvents} bubble events, "
                         f"source-off = {backEvents} bubble events")
            for label, n, b in empty:
                upperLimitRate = ZERO_COUNT_UPPER / (g["liveTimeSec"] / 60) if g["liveTimeSec"] > 0 else float("nan")
                lines.append(f"    mult {label:<3}: source-on = {n} events, source-off = {b} events"
                             f"  -> expected background {b * scale:0.3f} events,"
                             f" 68% upper limit {upperLimitRate:0.4f} count/min")
        if emptyBins:
            lines.append(f"  {emptyBins} of {allBins} bins have zero source-on events")
        else:
            lines.append(f"  none -- all {allBins} bins have at least one source-on event")
    report_section(GROUP_REPORT, "Multiplicity bins with zero source-on events", lines)


report_empty_bins([
    ("withoutDomeCut (pre-cut, Z > %scm, sim off)" % DOME_Z_THRESHOLD_CM_PRE, groupsWithoutDomeCut),
    ("withDomeCut (post-cut, Z > %scm, sim on)" % DOME_Z_THRESHOLD_CM_POST, groupsWithDomeCut),
])


# the source-activity normalization: the factor the *plotted* Seitz curve would have to be multiplied by
# for the data to line up with it. Everything here is in the units the plots are drawn in -- g["seitzRate"]
# is the blue curve, already scaled by build_groups()'s shared normalizationFactor, and the observed rates
# are the same background-subtracted numbers the points sit at -- so a ratio of 1.00 means "already
# aligned" and the multiplicity-1 column is the single-bubble number that sets the source activity.
#
# Rates are turned back into counts (rate * livetime) before pooling so groups combine by exposure rather
# than each short group counting as much as a long one. The observed side uses the unfloored rates: the
# drawn points are floored at 0, and averaging in a floored deficit bin would bias the ratio upward.
def report_source_activity_normalization(groupLists):
    lines = ["a ratio of 1.00 means already aligned. normalizationFactor here is FITTED (livetime-weighted",
             "per pipeline, over all 3 bins per group), not derived from the mult 1 column alone, so even",
             "the mult 1 ratio is not forced to exactly 1.00 -- how far it sits from 1.00 shows how much",
             "the fit's per-group total is shaped by the other multiplicity bins"]
    for pipelineLabel, groups, normalizationFactor in groupLists:
        lines += ["", f"=== {pipelineLabel} ==="]
        lines.append(f"    sim counts -> plotted rate factor currently in use "
                     f"(fitted, livetime-weighted): "
                     f"{normalizationFactor:0.6g} count/min per sim count")

        binIndices = range(len(rebinnedBinLabels))
        observedTotals = [sum(g["backSubUnfloored"][i] * g["liveTime"] for g in groups) for i in binIndices]
        plottedTotals = [sum(g["seitzRate"][i] * g["liveTime"] for g in groups) for i in binIndices]
        ratios = [o / p if p > 0 else float("nan") for o, p in zip(observedTotals, plottedTotals)]

        lines.append("  observed vs the plotted Seitz curve, pooled over every group [counts]:")
        lines.append(f"    {'':<12}" + "".join(f"{('mult ' + label):>12}" for label in rebinnedBinLabels))
        lines.append(f"    {'observed':<12}" + "".join(f"{v:12.2f}" for v in observedTotals))
        lines.append(f"    {'plotted sim':<12}" + "".join(f"{v:12.2f}" for v in plottedTotals))
        lines.append(f"    {'obs/sim':<12}" + "".join(f"{v:12.3f}" for v in ratios))

        singlesRatio = ratios[0]
        if np.isfinite(singlesRatio):
            lines.append(f"  source-activity normalization from the single-bubble count: {singlesRatio:0.4f}")
            rescaled = "   ".join(f"mult {label} -> {r / singlesRatio:0.3f}"
                                  for label, r in zip(rebinnedBinLabels[1:], ratios[1:]))
            lines.append(f"    scale the plotted sim curve by this and mult 1 lines up, leaving:  {rescaled}")

        # the same singles ratio group by group -- one shared activity number is only defensible if these
        # scatter no more than their own counting statistics
        lines.append("  per-group single-bubble obs/sim:")
        perGroup = []
        for g in groups:
            plotted = g["seitzRate"][0] * g["liveTime"]
            if plotted <= 0:
                continue
            ratio = g["backSubUnfloored"][0] * g["liveTime"] / plotted
            perGroup.append(ratio)
            lines.append(f"      Seitz {g['seitz']:6.2f} keV  ({group_pt_description(g)}):  {ratio:0.4f}   "
                         f"(livetime {g['liveTime']:0.2f} min, {g['binCountsRaw'][0]} source-on singles)")
        if len(perGroup) > 1:
            arr = np.asarray(perGroup)
            lines.append(f"    spread: min {arr.min():0.4f}, max {arr.max():0.4f}, mean {arr.mean():0.4f}, "
                         f"std {arr.std(ddof=1):0.4f} ({100 * arr.std(ddof=1) / arr.mean():0.1f}% of mean)")
    report_section(SUMMARY_REPORT,
                   "Source-activity normalization: the factor the plotted Seitz curve would have to be "
                   "multiplied by for the data to line up with it", lines)


report_source_activity_normalization([
    ("withoutDomeCut (pre-cut, Z > %scm, sim off)" % DOME_Z_THRESHOLD_CM_PRE,
     groupsWithoutDomeCut, normalizationFactorWithoutDomeCut),
    ("withDomeCut (post-cut, Z > %scm, sim on)" % DOME_Z_THRESHOLD_CM_POST,
     groupsWithDomeCut, normalizationFactorWithDomeCut),
])


# how far the plotted sim curve sits from the plotted data, over exactly the panels one figure shows.
# Both series are the ones actually drawn: g["backSub"] (floored at 0) and g["seitzRate"].
#
# Rates become counts (rate * livetime) before pooling, so a long group weighs more than a short one, and
# sigma puts the gap in units of the drawn error bar. normalizationFactor here IS fitted to the data
# (livetime-weighted per pipeline, over each group's own 3-bin total), so the "all bins" column is closer
# to a near-zero identity than neutronExcluded.py's -- but not exactly zero, since the fit matches each
# group's own total, not the grand total pooled across groups here, and not the per-bin shape at all.
def report_data_sim_percent_difference(groups, figureName):
    labels = rebinnedBinLabels
    pct = lambda data, sim: 100 * (data - sim) / sim if sim else float("nan")
    lines = [f"figure: {figureName}   ({len(groups)} groups, {len(groups) * len(labels)} bins)",
             "percent difference = 100 * (data - sim) / sim   [sim = the plotted Seitz curve]",
             "sigma = (data - sim) / upper error bar",
             "",
             "  per bin:",
             f"    {'Seitz':>7} {'T[min]':>8} " + "".join(f"{('mult ' + l):>18}" for l in labels)]
    for g in groups:
        cells = ""
        for data, sim, errHigh in zip(g["backSub"], g["seitzRate"], g["errHigh"]):
            sigma = (data - sim) / errHigh if errHigh else float("nan")
            cells += f"{pct(data, sim):+11.1f}% {sigma:+5.1f}s"
        lines.append(f"    {g['seitz']:7.2f} {g['liveTime']:8.2f} {cells}")

    dataTotals = [sum(g["backSub"][i] * g["liveTime"] for g in groups) for i in range(len(labels))]
    simTotals = [sum(g["seitzRate"][i] * g["liveTime"] for g in groups) for i in range(len(labels))]
    lines += ["",
              "  pooled over every panel [counts]:",
              f"    {'':<10}" + "".join(f"{('mult ' + l):>14}" for l in labels) + f"{'all bins':>14}",
              f"    {'data':<10}" + "".join(f"{v:14.2f}" for v in dataTotals) + f"{sum(dataTotals):14.2f}",
              f"    {'sim':<10}" + "".join(f"{v:14.2f}" for v in simTotals) + f"{sum(simTotals):14.2f}",
              f"    {'% diff':<10}" + "".join(f"{pct(d, m):+13.1f}%" for d, m in zip(dataTotals, simTotals))
              + f"{pct(sum(dataTotals), sum(simTotals)):+13.1f}%"]

    perBin = [pct(data, sim) for g in groups for data, sim in zip(g["backSub"], g["seitzRate"])]
    arr = np.asarray([v for v in perBin if np.isfinite(v)])
    if len(arr):
        lines.append(f"    spread over the {len(arr)} bins: mean {arr.mean():+0.1f}%, "
                     f"median {np.median(arr):+0.1f}%, range {arr.min():+0.1f}% to {arr.max():+0.1f}%, "
                     f"mean |diff| {np.abs(arr).mean():0.1f}%")

    # a column that lands on the same side in every panel is a systematic, not scatter -- worth saying so
    # explicitly, since no single bin need reach 2 sigma for the pattern across panels to be real
    for i, label in enumerate(labels):
        column = [pct(g["backSub"][i], g["seitzRate"][i]) for g in groups]
        if all(v > 0 for v in column) or all(v < 0 for v in column):
            direction = "above" if column[0] > 0 else "below"
            lines.append(f"    note: mult {label} is {direction} the sim in all {len(groups)} panels "
                         f"({min(column):+0.1f}% to {max(column):+0.1f}%) -- coherent, not scatter")
    report_section(SUMMARY_REPORT, "Data vs plotted sim, percent difference", lines)


report_data_sim_percent_difference(groupsWithDomeCut[FINAL_PLOT_SLICE], "combinedMultiplicityFinal.pdf")


# the rate integrated over every multiplicity bin, one line per Seitz threshold -- what a group's total
# bubble rate is, without the multiplicity structure the other reports slice it by.
#
# The data error is re-derived from the SUMMED raw counts rather than quadrature-summing the per-bin
# errors. Same reasoning _rate_fields() applies when it re-derives the rebinned bins from raw counts: the
# zero-count upper limit is a property of an empty bin, so adding several of them in quadrature would give
# a well-populated total a much larger error than its own count implies. Unfloored, like the fits use.
def report_integrated_rates(groupLists):
    lines = ["integrated = summed over every multiplicity bin (1..5+), count/min",
             "data errors re-derived from the summed raw counts, not quadrature-summed across bins",
             "sim = the plotted Seitz curve, summed the same way"]
    for pipelineLabel, groups in groupLists:
        lines += ["", f"=== {pipelineLabel} ===",
                  f"    {'Seitz':>7} {'T[min]':>8}   {'data [count/min]':>26} {'sim':>9} {'data/sim':>9}"
                  f"   {'source-on':>9}  setpoints"]
        for g in groups:
            _, _, _, subBins, errLow, errHigh, _, _ = background_subtract(
                [sum(g["binCountsRaw"])], g["liveTimeSec"],
                [sum(g["backCountsRaw"])], g["backgroundTimeSec"])
            liveMin = g["liveTime"]
            data, low, high = subBins[0] / liveMin, errLow[0] / liveMin, errHigh[0] / liveMin
            sim = sum(g["seitzRate"])
            ratio = data / sim if sim else float("nan")
            lines.append(f"    {g['seitz']:7.2f} {liveMin:8.2f}   "
                         f"{data:9.4f} (+{high:0.4f}/-{low:0.4f}) {sim:9.4f} {ratio:9.3f}"
                         f"   {sum(g['binCountsRaw']):9d}  {group_pt_description(g)}")
    report_section(SUMMARY_REPORT, "Integrated rate per Seitz threshold", lines)


report_integrated_rates([
    ("withoutDomeCut (pre-cut, Z > %scm, sim off)" % DOME_Z_THRESHOLD_CM_PRE, groupsWithoutDomeCut),
    ("withDomeCut (post-cut, Z > %scm, sim on)" % DOME_Z_THRESHOLD_CM_POST, groupsWithDomeCut),
])

## plot making
"""
COMBINED PAPER PLOT -- side by side, dome cut vs no dome cut
"""
# groupsWithCut/groupsWithoutCut must be the same length, in the same (p,T) order (guaranteed since both come from build_groups() over the same pToUse, sorted by seitz)
def plot_combined_multiplicity_comparison(groupsWithCut, groupsWithoutCut, savepath, groupsPerRow=3):
    binLabels = rebinnedBinLabels
    numBins = len(binLabels)
    barWidth = 0.9
    pairGap = 0.5   # gap between the "cut" and "no cut" halves of one threshold
    gap = 1.0        # gap between different thresholds
    colWidthInches = 5.0

    nGroups = len(groupsWithCut)
    nRows = int(np.ceil(nGroups / groupsPerRow))
    fig, axes = plt.subplots(nRows, 1, figsize=(colWidthInches, 3.4 * nRows), squeeze=False)
    axes = axes[:, 0]

    globalMax = max(
        max(max(g["seitzRate"]), max(b + e for b, e in zip(g["backSub"], g["errHigh"])))
        for groupList in (groupsWithCut, groupsWithoutCut) for g in groupList
    )

    for rowIdx, ax in enumerate(axes):
        trans = ax.get_xaxis_transform()
        rowWith = groupsWithCut[rowIdx * groupsPerRow:(rowIdx + 1) * groupsPerRow]
        rowWithout = groupsWithoutCut[rowIdx * groupsPerRow:(rowIdx + 1) * groupsPerRow]

        pos = 0
        for gi, (gWith, gWithout) in enumerate(zip(rowWith, rowWithout)):
            xsWith = np.arange(pos, pos + numBins)
            pos += numBins + pairGap
            xsWithout = np.arange(pos, pos + numBins)
            pos += numBins

            ax.bar(xsWith, gWith["seitzRate"], width=barWidth, color="lightblue", edgecolor="steelblue", zorder=1)
            ax.errorbar(xsWith, gWith["backSub"], yerr=[gWith["errLow"], gWith["errHigh"]], fmt='o', color="red",
                        ecolor="red", zorder=2, markersize=3, elinewidth=1, capsize=2)

            ax.bar(xsWithout, gWithout["seitzRate"], width=barWidth, color="lightblue", edgecolor="steelblue",
                   hatch="//", zorder=1)
            ax.errorbar(xsWithout, gWithout["backSub"], yerr=[gWithout["errLow"], gWithout["errHigh"]], fmt='^',
                        color="darkorange", ecolor="darkorange", zorder=2, markersize=3, elinewidth=1, capsize=2)

            for x, label in zip(xsWith, binLabels):
                ax.text(x, -0.05, label, transform=trans, ha='center', va='top', fontsize=9)
            for x, label in zip(xsWithout, binLabels):
                ax.text(x, -0.05, label, transform=trans, ha='center', va='top', fontsize=9)
            ax.text((xsWith[0] + xsWith[-1]) / 2, -0.14, "cut", transform=trans, ha='center', va='top',
                    fontsize=8, style='italic')
            ax.text((xsWithout[0] + xsWithout[-1]) / 2, -0.14, "no cut", transform=trans, ha='center', va='top',
                    fontsize=8, style='italic')

            # seitz threshold label (same for both halves -- seitz doesn't depend on the dome cut)
            center = (xsWith[0] + xsWithout[-1]) / 2
            ax.text(center, 0.97, f'{gWith["seitz"]:0.2f}', transform=trans, ha='center', va='top', fontsize=10)

            # seperator
            if gi < len(rowWith) - 1:
                ax.axvline(pos + gap / 2 - 0.5, linestyle='--', linewidth=0.7, color='gray', zorder=0)

            pos += gap

        ax.set_ylim(0, globalMax * 1.25)
        ax.set_xlim(-1, pos - gap)
        ax.set_xticks([])
        ax.tick_params(axis='y', labelsize=12)

    axes[0].set_title(r'$Q_{seitz}$ [keV]', loc='left', fontsize=16, pad=2)
    axes[0].set_title(r'$^{252}$Cf Source: dome cut vs no cut', loc='right', fontsize=12, pad=2)
    axes[nRows // 2].set_ylabel("Rate [count/min]", fontsize=16)
    axes[-1].set_xlabel("Bubble Multiplicity", fontsize=12, labelpad=40)

    legendHandles = [
        plt.Line2D([0], [0], marker='o', linestyle='', color='red', label='Data (dome cut)'),
        plt.Line2D([0], [0], marker='^', linestyle='', color='darkorange', label='Data (no cut)'),
        plt.Rectangle((0, 0), 1, 1, facecolor='lightblue', edgecolor='steelblue', label='Seitz pred. (dome cut)'),
        plt.Rectangle((0, 0), 1, 1, facecolor='lightblue', edgecolor='steelblue', hatch='//', label='Seitz pred. (no cut)'),
    ]
    fig.legend(handles=legendHandles, loc='upper center', ncol=2, fontsize=9, bbox_to_anchor=(0.5, 1.0))

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.3, top=0.88)
    fig.savefig(savepath)
    plt.close(fig)

# same grid-of-thresholds layout as plot_combined_multiplicity_comparison() above, but a single
# series (no cut/no-cut pairing) -- pulls the standalone data for just one side out on its own
def plot_combined_multiplicity_single(groups, savepath, negLnLByGroup, chi2ByGroup, sourceLabel,
                                       groupsPerRow=4, showPerGroupStats=False):
    binLabels = rebinnedBinLabels
    numBins = len(binLabels)
    barWidth = 1.0
    gap = 0.6
    colWidthInches = 3.5

    nRows = int(np.ceil(len(groups) / groupsPerRow))
    fig, axes = plt.subplots(nRows, 1, figsize=(colWidthInches, 3.2 * nRows), squeeze=False)
    axes = axes[:, 0]

    globalMax = max(
        max(max(g["seitzRate"]), max(b + e for b, e in zip(g["backSub"], g["errHigh"])))
        for g in groups
    )

    for rowIdx, ax in enumerate(axes):
        trans = ax.get_xaxis_transform()
        rowGroups = groups[rowIdx * groupsPerRow:(rowIdx + 1) * groupsPerRow]

        rowIdx0 = rowIdx * groupsPerRow
        pos = 0
        for gi, g in enumerate(rowGroups):
            xs = np.arange(pos, pos + numBins)

            ax.bar(xs, g["seitzRate"], width=barWidth, color="lightblue", edgecolor="steelblue", zorder=1)
            ax.errorbar(xs, g["backSub"], yerr=[g["errLow"], g["errHigh"]], fmt='o', color="red",
                        ecolor="red", zorder=2, markersize=3, elinewidth=1, capsize=2)

            for x, label in zip(xs, binLabels):
                ax.text(x, -0.05, label, transform=trans, ha='center', va='top', fontsize=12)

            center = (xs[0] + xs[-1]) / 2
            ax.text(center, 0.97, f'{g["seitz"]:0.2f}', transform=trans, ha='center', va='top', fontsize=10)

            if showPerGroupStats:
                groupNegLnL = negLnLByGroup[rowIdx0 + gi]
                groupChi2 = chi2ByGroup[rowIdx0 + gi]
                perGroupText = r'$-2\Delta\ln\mathcal{L}$' + f'={groupNegLnL:0.1f}'
                ax.text(center, 0.86, perGroupText, transform=trans, ha='center', va='top', fontsize=6)

            if gi < len(rowGroups) - 1:
                ax.axvline(pos + numBins + gap / 2 - 0.5, linestyle='--', linewidth=0.7,
                           color='gray', zorder=0)

            pos += numBins + gap

        ax.set_ylim(0, globalMax * 1.2)
        ax.set_xlim(-1, pos - gap)
        ax.set_xticks([])
        ax.tick_params(axis='y', labelsize=12)

    # the combined totals aren't drawn on the figure anymore -- written to captionStats.txt so they can go in the caption
    totalNegLnL = sum(negLnLByGroup)
    totalChi2 = sum(chi2ByGroup)
    totalBins = len(groups) * numBins
    report_section(SUMMARY_REPORT, "Figure fit statistics",
                   [f"{os.path.basename(savepath)}: "
                    f"-2*deltaLnL/n.d.o.f. = {totalNegLnL:0.1f}/({totalBins} - 1), "
                    f"chi2/n.d.o.f. = {totalChi2:0.1f}/({totalBins} - 1)"])

    axes[0].set_title(r'$Q_{seitz}$ [keV]', loc='left', fontsize=16, pad=2)
    axes[0].set_title(sourceLabel, loc='right', fontsize=12, pad=2)
    # figure-level label so it sits centered on the whole stack of rows rather than on whichever
    # single row happens to be the middle one (which is off-center whenever nRows is even)
    fig.supylabel("Rate [count/min]", fontsize=16)
    axes[-1].set_xlabel("Bubble Multiplicity", fontsize=12, labelpad=28)

    legendHandles = [
        plt.Line2D([0], [0], marker='o', linestyle='', color='red', label='Data (Source ON - OFF)'),
        plt.Rectangle((0, 0), 1, 1, facecolor='lightblue', edgecolor='steelblue', label='Source MC'),
    ]
    fig.legend(handles=legendHandles, loc='upper center', ncol=2, fontsize=9, bbox_to_anchor=(0.5, 1.0))

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.15, top=0.9 if nRows == 1 else 0.92)

    fig.savefig(savepath, bbox_inches="tight")
    plt.close(fig)

plot_combined_multiplicity_comparison(
    groupsWithDomeCut, groupsWithoutDomeCut,
    savepath=output_path("comparison", "combinedMultiplicityComparison.png"),
)

# same plot, but the "no cut" side drops the PRE baseline too -- dome cut vs truly nothing, for
# showing how much the cuts (PRE included) actually move the data
plot_combined_multiplicity_comparison(
    groupsWithDomeCut, groupsNoCutAtAll,
    savepath=output_path("comparison", "combinedMultiplicityComparisonNoCutAtAll.png"),
)

# standalone (single-series) version of each side that appears in the two comparison plots above --
# calls moved below fit_threshold_series() (search "STANDALONE COMBINED MULTIPLICITY PLOTS"), since they
# need neg_ln_l_calc()/global_normalization_factor(), which aren't defined yet at this point in the script

"""
1-BUBBLE Z DISTRIBUTION, FOUR LOWEST SEITZ THRESHOLDS
"""
# how many of the lowest-seitz groups to plot
numLowestSeitzZDist = 4

# {(p, T): [z, z, ...]} -- real reco Z for confirmed single-bubble events, grouped by (p, T); drops dome events unless applyDomeExclusion=False (then kept, only counted for the print below)
def load_single_bubble_z_by_group(runList, applyDomeExclusion=True):
    byGroup = {}
    recoLookupCache = {}
    evDataCache = {}
    regionExcludedCount = 0
    zExcludedCount = 0
    for run, ev, region in iter_single_bubble_events(HANDSCAN_DIR, runList):
        if run not in recoLookupCache:
            recoPath = os.path.join(RECON_DIR, run, "reco.sbc")
            recoLookupCache[run] = build_reco_lookup(Streamer(recoPath).to_dict()) if os.path.exists(recoPath) else None
        recoLookup = recoLookupCache[run]
        if recoLookup is None:
            continue
        coord = first_valid_coord(recoLookup, ev)
        if coord is None:
            continue
        if "dome" in excludedRegions and is_real_dome_event(region, coord):
            if region == 2:
                regionExcludedCount += 1
            else:
                zExcludedCount += 1
            if applyDomeExclusion:
                continue
        pset = _event_pset_temp(run, ev, evDataCache)
        if pset is None:
            continue
        lo, hi, temp = pset
        if lo != hi:
            continue
        byGroup.setdefault((lo, temp), []).append(float(coord[2]))
    return byGroup


# z values below this, and r^2 values above this, are non-physical, so leave them out of the distribution plots
MIN_PHYSICAL_Z_MM = -150 * 10
MAX_PHYSICAL_R2_MM2 = 20000


# source/background counts per Z-bin, each divided by its own livetime [minutes] -> count/min
def plot_z_distribution(sourceZ, backgroundZ, seitz, savepath, sourceLiveTime, backgroundLiveTime, numBins=50 // 4):
    sourceZ = [z for z in sourceZ if z >= MIN_PHYSICAL_Z_MM]
    backgroundZ = [z for z in backgroundZ if z >= MIN_PHYSICAL_Z_MM]
    allZ = sourceZ + backgroundZ
    edges = np.linspace(min(allZ), max(allZ), numBins + 1) if allZ else numBins

    backgroundCounts, _ = np.histogram(backgroundZ, bins=edges)
    sourceCounts, _ = np.histogram(sourceZ, bins=edges)
    backgroundRate = backgroundCounts / backgroundLiveTime if backgroundLiveTime > 0 else np.zeros_like(backgroundCounts, dtype=float)
    sourceRate = sourceCounts / sourceLiveTime

    plt.figure(figsize=(8, 6))
    plt.stairs(backgroundRate, edges, color="gray", linewidth=1.5, label="Background")
    plt.stairs(sourceRate, edges, color="steelblue", linewidth=1.5, label="Source")
    plt.axvline(DOME_Z_THRESHOLD_CM * 10, color="red", linestyle="--",
                label=f"dome cut ({DOME_Z_THRESHOLD_CM * 10:.0f} mm)")
    plt.xlabel("Z [mm]", fontsize=16)
    plt.ylabel("Rate [count/min]", fontsize=16)
    plt.title(f"1-bubble Z distribution, Seitz = {seitz:0.2f} keV", fontsize=14)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig(savepath)
    plt.close()


singleBubbleZByGroup = load_single_bubble_z_by_group(neutronRuns, applyDomeExclusion=False)
# background Z split the same way the subtraction is, so each group's overlay is background from its own (p, T)
backgroundZByGroup = load_single_bubble_z_by_group(backgroundList, applyDomeExclusion=False)

# detector wall + dome boundary, r^2-vs-z view -- same construction as draw_r2z_guides() in
# ../SBC_handscan/reconAna/reconAna.py, inlined here (rather than imported) since that script
# lives outside this project and isn't guaranteed to be present wherever this pipeline runs
def draw_detector_r2z_guides(ax):
    ax.vlines((25.4 * 4.525) ** 2, 25.4 * -8.75, 25.4 * (14.71997 - 15.358), color='r')
    ax.vlines((25.4 * 4.725) ** 2, ymin=25.4 * -8.75, ymax=25.4 * (14.71997 - 15.358), color='r')

    theta = np.linspace(0, 1.19367, 400)
    rcirc = 2 * 25.4
    ax.plot((rcirc * np.cos(theta) + 25.4 * 2.725) ** 2,
            rcirc * np.sin(theta) + 25.4 * (14.71997 - 15.358), c='r')

    rcirc = 1.8 * 25.4
    ax.plot((rcirc * np.cos(theta) + 25.4 * 2.725) ** 2,
            rcirc * np.sin(theta) + 25.4 * (14.71997 - 15.358), c='r')


# scatter of every confirmed single-bubble event with a valid 3D reco coord: Z vs r^2 = X^2+Y^2 (mm^2),
# with the pre-cut dome Z threshold and the detector's own wall/dome boundary overlaid --
# unfiltered (applyDomeExclusion=False) so the plot shows where the cut line actually falls
# relative to the full reconstructed event cloud
def plot_z_vs_r2_distribution(sourcePositions, backgroundPositions, savepath):
    def r2_and_z(positions):
        physical = [
            pos for pos in positions
            if pos[4] >= MIN_PHYSICAL_Z_MM and pos[2] ** 2 + pos[3] ** 2 <= MAX_PHYSICAL_R2_MM2
        ]
        return [pos[2] ** 2 + pos[3] ** 2 for pos in physical], [pos[4] for pos in physical]

    backgroundR2, backgroundZ = r2_and_z(backgroundPositions)
    sourceR2, sourceZ = r2_and_z(sourcePositions)

    plt.figure(figsize=(8, 6))
    plt.scatter(backgroundR2, backgroundZ, s=4, color="gray", alpha=0.4, label="Background")
    plt.scatter(sourceR2, sourceZ, s=4, color="steelblue", alpha=0.4, label=r"$^{252}$Cf")
    plt.axhline(DOME_Z_THRESHOLD_CM_PRE * 10, color="darkorange", linestyle="--")
    draw_detector_r2z_guides(plt.gca())
    plt.xlabel(r"$r^2$ [mm$^2$]", fontsize=16)
    plt.ylabel("Z [mm]", fontsize=16)
    plt.title("1-bubble events with 3D reco: Z vs $r^2$", fontsize=14)
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.savefig(savepath)
    plt.close()

plot_z_vs_r2_distribution(
    load_single_bubble_positions(neutronRuns, label="source (Z-vs-r^2)", applyDomeExclusion=False),
    load_single_bubble_positions(backgroundList, label="background (Z-vs-r^2)", applyDomeExclusion=False),
    savepath=output_path("comparison", "zVsR2Distribution.png"),
)

# if true, use 1,2,3+ if false use 1,2,3,4,5+
useRebinnedThresholdPlots = True

"""
THEORETICAL THRESHOLDS
"""
# range to check for ratio matching
thresholdRange = np.arange(0.2, 30.1, 0.1)

# range of multipliers to scan for the single-A "threshold = A * seitz" fit below
seitzMultiplierRange = np.arange(0.1, 3.01, 0.01)

# fit normalization mode: if True, scale every group's predicted counts by one shared factor across the whole dataset; if False, scale each to its own observed total (like oldCode.py)
useGlobalNormalization = True

# if True, also compute the *other* normalization mode and generate the *Compare* plots overlaying both. If False (default), just use useGlobalNormalization above
compareNormalizationModes = False

# if True, also fit/plot the no-singles (2, 3+ only) variant alongside the normal fit. Off by default since it's almost never needed
computeNoSinglesFit = False

# master switch for every threshold-fitting step: the per-group best-fit scan over thresholdRange, the
# avg-group fit, the no-singles fit, and the single-A "threshold = A * seitz" fit, plus every plot that
# needs a fitted threshold (negLnLScans/, fitToSeitzRatio*, linhistFit*, avgseitzFit*, linHistsA/).
# Setting this False skips all of that (the slow part) but leaves everything else intact: rates,
# Z distributions, plain linhists, and the combined multiplicity plots with their negLnL/chi2 stats
# evaluated at each group's own Seitz threshold
computeBestFit = False

def group_observed_total(dataGroup):
    rate = dataGroup["backSubUnfloored"] if useRebinnedThresholdPlots else dataGroup["backSubUnflooredFull"]
    return sum(v * dataGroup["liveTime"] for v in rate)

# livetime-weighted average of each group's own (sim-predicted total / observed total) ratio -- same
# weighting pattern build_groups() uses for normalizationFactor, but here expressed relative to each
# group's own observed total so it can be applied as targetTotal = globalNormFactor * sum(observed) in
# neg_ln_l_calc()/chi_squared_diagnostic(). g["seitzRate"] is already normalizationFactor-calibrated, so
# this comes out close to 1 (deviating only by genuine group-to-group scatter) instead of always >= 1 --
# the old version compared each group's observed total only to the *other groups'* observed totals
# (sum(w^2)/sum(w)/mean(w), which is >= 1 by Cauchy-Schwarz whenever livetimes/rates differ across
# groups) and never referenced the simulation at all, so it silently inflated every predicted total
# NOTE: this is the last livetime-weighted, data-fitted quantity in this variant. It is inert as shipped --
# computeBestFit is False, and at a group's own Seitz Q neg_ln_l_calc()/chi_squared_diagnostic() bypass it
# and read g["seitzRate"] directly -- so it touches no output. Turning computeBestFit on WOULD reintroduce
# both livetime weighting and normalization-to-data into the threshold scans. Left as-is rather than
# silently changed, so the fits still behave the way neutronExcluded.py's do if you enable them
def global_normalization_factor(dataGroups):
    weights = [g["liveTime"] for g in dataGroups]
    ratioByGroup = [
        sum(g["seitzRate"]) * g["liveTime"] / group_observed_total(g) for g in dataGroups
    ]
    return sum(w * r for w, r in zip(weights, ratioByGroup)) / sum(weights)

def normalization_mode_label(useGlobalNorm):
    return "Global normalization" if useGlobalNorm else "Per-threshold normalization"

# Poisson profile-likelihood test statistic for one bin -- replaces (obs-pred)^2/err^2.
# -2*ln[L(S, mu_b_hat)/L_saturated], profiling out the background rate mu_b_hat(S) instead
# of propagating obs/pred through a hand-built asymmetric Gaussian error. N = raw source-run
# count, B = raw background-run count, S = predicted signal count, tau = backgroundTime /
# sourceLiveTime (same units, e.g. both seconds). Standard on/off counting-experiment profile
# likelihood, e.g. Cowan, Cranmer, Gross & Vitells (2011), "Asymptotic formulae...".
def neg2_delta_ln_l(N, B, S, tau):
    onePlusTau = 1 + tau
    discriminant = (onePlusTau * S - N - B) ** 2 + 4 * onePlusTau * B * S
    muB = ((N + B) - onePlusTau * S + np.sqrt(discriminant)) / (2 * onePlusTau)
    muB = max(muB, 0.0)  # guards float noise right at the B=0/S=0 boundary

    sPlusMuB = S + muB
    # x*ln(x) -> 0 as x -> 0, by convention (avoids 0*log(0) for empty bins)
    nTerm = N * np.log(N / sPlusMuB) if N > 0 and sPlusMuB > 0 else 0.0
    bTerm = B * np.log(B / (tau * muB)) if B > 0 and tau * muB > 0 else 0.0

    return 2 * (nTerm + bTerm + sPlusMuB - N + tau * muB - B)

# globalNormFactor and simExcludedRegions are per-side, since they depend on which groups list / sim cut is active.
# Named generically (neg_ln_l_calc, not neg2_delta_ln_l_calc) because the rest of the fit/plot pipeline --
# compute_best_fit, neg_ln_l_confidence_interval, plot_neg_ln_l_scan, ... -- just minimizes whatever curve
# this returns; the actual -2*deltaLnL(S) math lives in neg2_delta_ln_l() above.
def neg_ln_l_calc(dataGroup, estThreshold, globalNormFactor, simExcludedRegions, useGlobalNorm=None):
    if useGlobalNorm is None:
        useGlobalNorm = useGlobalNormalization

    if useRebinnedThresholdPlots:
        rate = dataGroup["backSubUnfloored"]
        predictedCounts = rebin(get_multiplicity_counts(estThreshold, excludedRegionsOverride=simExcludedRegions)[0])
        N = rebin(dataGroup["binCountsRaw"])
        B = rebin(dataGroup["backCountsRaw"])
    else:
        rate = dataGroup["backSubUnflooredFull"]
        predictedCounts = get_multiplicity_counts(estThreshold, excludedRegionsOverride=simExcludedRegions)[0]
        N = dataGroup["binCountsRaw"]
        B = dataGroup["backCountsRaw"]

    liveTime = dataGroup["liveTime"]

    # at a group's own Seitz threshold (how every combinedMultiplicity*/PerGroupStats plot's
    # negLnLByGroup calls this), use dataGroup["seitzRate"] directly -- same plotted curve
    # chi_squared_diagnostic() now matches -- converted to counts (*liveTime) since N/B/S must share
    # units for neg2_delta_ln_l(). Away from that exact threshold (fit_groups()'s thresholdRange scan,
    # run_avg_group_pipeline()'s avgGroup which has no "seitzRate") there's no single plotted curve to
    # match, so fall back to the shape-normalized-to-targetTotal curve used to actually drive the fit
    if useRebinnedThresholdPlots and "seitzRate" in dataGroup and estThreshold == dataGroup.get("seitz"):
        predicted = [r * liveTime for r in dataGroup["seitzRate"]]
    else:
        # target total is still set from the background-subtracted rate -- only the goodness-of-fit
        # comparison below switched to raw counts, so the predicted curve's overall scale is unchanged
        observed = [v * liveTime for v in rate]
        targetTotal = globalNormFactor * sum(observed) if useGlobalNorm else sum(observed)
        predicted = seitz_count(counts_to_ratios(predictedCounts), targetTotal)

    tau = dataGroup["backgroundTimeSec"] / dataGroup["liveTimeSec"]
    return sum(neg2_delta_ln_l(n, b, s, tau) for n, b, s in zip(N, B, predicted))

# chi-squared analog of neg_ln_l_calc() above -- the (obs-pred)^2/err^2 sum this pipeline used before
# switching to the Poisson profile likelihood, kept only for the neg_ln_l-vs-chi2 diagnostic plot.
# Same shape/target-total logic as neg_ln_l_calc(), just compared via a hand-built asymmetric Gaussian
# error instead of neg2_delta_ln_l()'s exact Poisson treatment
def chi_squared_diagnostic(dataGroup, estThreshold, globalNormFactor, simExcludedRegions, useGlobalNorm=None,
                            debug=False):
    if useGlobalNorm is None:
        useGlobalNorm = useGlobalNormalization

    if useRebinnedThresholdPlots:
        binLabels = rebinnedBinLabels
        rate, errLowRate, errHighRate = (dataGroup["backSubUnfloored"], dataGroup["errLowUnfloored"],
                                         dataGroup["errHighUnfloored"])
        predictedCounts = rebin(get_multiplicity_counts(estThreshold, excludedRegionsOverride=simExcludedRegions)[0])
    else:
        binLabels = ["1", "2", "3", "4", "5+"]
        rate = dataGroup["backSubUnflooredFull"]
        errLowRate, errHighRate = dataGroup["errLowUnflooredFull"], dataGroup["errHighUnflooredFull"]
        predictedCounts = get_multiplicity_counts(estThreshold, excludedRegionsOverride=simExcludedRegions)[0]

    # at a group's own Seitz threshold (how every combinedMultiplicity*/PerGroupStats plot calls this),
    # use dataGroup["seitzRate"] directly -- the exact same curve drawn as the plot's blue bars -- instead
    # of re-deriving a separately-normalized prediction via globalNormFactor. Away from that threshold
    # (the A-multiplier scan in chi_squared_diagnostic_seitz_multiplier, A != 1) there's no plotted curve
    # to match at that threshold, so fall back to the shape-normalized-to-targetTotal curve
    if useRebinnedThresholdPlots and "seitzRate" in dataGroup and estThreshold == dataGroup.get("seitz"):
        predictedRate = dataGroup["seitzRate"]
        predictedRateSource = "seitzRate (plotted curve)"
    else:
        targetTotalRate = globalNormFactor * sum(rate) if useGlobalNorm else sum(rate)
        predictedRate = seitz_count(counts_to_ratios(predictedCounts), targetTotalRate)
        predictedRateSource = f"re-normalized, targetTotalRate={targetTotalRate:0.4f}"

    if debug:
        pT = group_debug_label(dataGroup)
        normalizedPredictedRate = [f'{p:0.4f}' for p in predictedRate]
        chi2_debug(f"  [chi2 debug] threshold={estThreshold:0.2f} keV  {pT}  liveTime={dataGroup['liveTime']:0.2f} min  "
                   f"predictedRate[{predictedRateSource}]={normalizedPredictedRate}")
        # raw event counts behind the subtraction, binned the same way as the chi2 terms below
        toDebugBins = rebin if useRebinnedThresholdPlots else list
        sourceOn = toDebugBins(dataGroup["binCountsRaw"])
        sourceOff = toDebugBins(dataGroup["backCountsRaw"])
        scale = dataGroup["backgroundScale"]
        chi2_debug(f"    raw counts: sourceOn={sourceOn}  sourceOff={sourceOff}  "
                   f"scaledSourceOff={[f'{c * scale:0.2f}' for c in sourceOff]}")
        plotted = dataGroup["backSub"] if useRebinnedThresholdPlots else dataGroup["backSubFull"]
        chi2_debug(f"    obsRate below is the unfloored subtraction (what the fit sees); plotted (floored at 0) = "
                   f"{[f'{v:0.4f}' for v in plotted]}")
        chi2_debug(f"    scale = sourceLiveTime/backgroundLiveTime = {dataGroup['liveTimeSec']:0.1f}s"
                   f"/{dataGroup['backgroundTimeSec']:0.1f}s = {scale:0.4f}")
        if dataGroup["backgroundTimeSec"] <= 0:
            chi2_debug("    WARNING: no background runs at this (p, T) -- nothing subtracted, and the fit's background rate is unconstrained here")

    chi2 = 0.0
    for label, obs, eLow, eHigh, pred in zip(binLabels, rate, errLowRate, errHighRate, predictedRate):
        err = eLow if obs >= pred else eHigh
        contribution = 0.0 if err == 0 else ((obs - pred) / err) ** 2
        chi2 += contribution
        if debug:
            chi2_debug(f"    mult {label:<3}: obsRate={obs:0.4f}  predRate={pred:0.4f}  "
                       f"errLow={eLow:0.4f}  errHigh={eHigh:0.4f}  errUsed={err:0.4f}  chi2contrib={contribution:0.4f}")
    if debug:
        chi2_debug(f"    -> chi2 total = {chi2:0.4f}")
    return chi2

# same idea as neg_ln_l_calc, but drops the multiplicity==1 bin and fits only the multi-bubble ones (2 and 3+, or just 2+) -- always rebinned, results go in g["bestFitNoSingles"], never g["bestFit"]
def neg_ln_l_calc_no_singles(dataGroup, estThreshold, simExcludedRegions):
    rate = dataGroup["backSubUnfloored"]
    predictedCounts = rebin(get_multiplicity_counts(estThreshold, excludedRegionsOverride=simExcludedRegions)[0])[1:]
    N = rebin(dataGroup["binCountsRaw"])[1:]
    B = rebin(dataGroup["backCountsRaw"])[1:]

    liveTime = dataGroup["liveTime"]
    observed = [v * liveTime for v in rate][1:]
    targetTotal = sum(observed)
    predicted = seitz_count(counts_to_ratios(predictedCounts), targetTotal)

    tau = dataGroup["backgroundTimeSec"] / dataGroup["liveTimeSec"]
    return sum(neg2_delta_ln_l(n, b, s, tau) for n, b, s in zip(N, B, predicted))


def neg_ln_l_confidence_interval(gridKev, negLnLCurve, bestIdx):
    """threshold values where the fit statistic crosses its minimum + 1, on each side of the best fit"""
    minNegLnL = negLnLCurve[bestIdx]
    target = minNegLnL + 1.0

    lowThreshold = gridKev[0]
    for i in range(bestIdx, 0, -1):
        if negLnLCurve[i - 1] >= target:
            x0, x1 = gridKev[i - 1], gridKev[i]
            y0, y1 = negLnLCurve[i - 1], negLnLCurve[i]
            frac = (target - y0) / (y1 - y0) if y1 != y0 else 0.0
            lowThreshold = x0 + frac * (x1 - x0)
            break

    highThreshold = gridKev[-1]
    for i in range(bestIdx, len(negLnLCurve) - 1):
        if negLnLCurve[i + 1] >= target:
            x0, x1 = gridKev[i], gridKev[i + 1]
            y0, y1 = negLnLCurve[i], negLnLCurve[i + 1]
            frac = (target - y0) / (y1 - y0) if y1 != y0 else 0.0
            highThreshold = x0 + frac * (x1 - x0)
            break

    return lowThreshold, highThreshold


# given a negLnL value per entry of thresholdRange, locates the best-fit threshold, its 1-sigma bounds, and the best-fit rate curve -- shared by every fit in the pipeline
# builds the standard {threshold, thresholdErrLow, thresholdErrHigh, negLnL, rate, rateFull} fit dict for an already-known threshold/negLnL -- rate is scaled to match dataGroup's own observed total, same as the 0keV/avg-seitz reference lines. Shared by compute_best_fit() below and fit_seitz_multiplier()'s per-group bestFitA
def build_fit_result(dataGroup, threshold, negLnL, simExcludedRegions, thresholdErrLow=0.0, thresholdErrHigh=0.0):
    bestFitRatios = counts_to_ratios(
        get_multiplicity_counts(threshold, excludedRegionsOverride=simExcludedRegions)[0]
    )
    bestFitRateFull = seitz_count(bestFitRatios, sum(dataGroup["backSubUnflooredFull"]))
    return {
        "threshold": threshold,
        "thresholdErrLow": thresholdErrLow,
        "thresholdErrHigh": thresholdErrHigh,
        "negLnL": negLnL,
        "rate": rebin(bestFitRateFull),
        "rateFull": bestFitRateFull,
    }


def compute_best_fit(dataGroup, negLnLCurve, simExcludedRegions):
    bestIdx = int(np.argmin(negLnLCurve))
    bestThreshold = thresholdRange[bestIdx]
    lowThreshold, highThreshold = neg_ln_l_confidence_interval(thresholdRange, negLnLCurve, bestIdx)

    fit = build_fit_result(
        dataGroup, bestThreshold, negLnLCurve[bestIdx], simExcludedRegions,
        thresholdErrLow=bestThreshold - lowThreshold, thresholdErrHigh=highThreshold - bestThreshold,
    )
    fit["thresholdRange"] = thresholdRange
    fit["negLnLCurve"] = negLnLCurve
    return fit

# -2*deltaLnL vs threshold scan
def plot_neg_ln_l_scan(gridKev, negLnLCurve, bestThreshold, bestNegLnL, savepath, xlabel="Threshold [keV]",
                        bestLabel=None, sigmaRange=None, infoText=None):
    if bestLabel is None:
        bestLabel = f"Best fit: {bestThreshold:0.2f} keV (" + r"$-2\Delta\ln L$" + f"={bestNegLnL:0.2f})"
    plt.figure(figsize=(8, 6))
    plt.plot(gridKev, negLnLCurve, 'o', markersize=3, color="steelblue")
    plt.axvline(bestThreshold, color='red', linestyle='--', label=bestLabel)
    if sigmaRange is not None:
        lowX, highX = sigmaRange
        plt.axvline(lowX, color='gray', linestyle=':', label=r"1$\sigma$ range: " + f"[{lowX:0.3f}, {highX:0.3f}]")
        plt.axvline(highX, color='gray', linestyle=':')
    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel(r"$-2\Delta\ln L$", fontsize=16)
    if infoText:
        plt.gca().text(0.02, 0.95, infoText, transform=plt.gca().transAxes, fontsize=11, va='top')
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig(savepath)
    plt.close()

# generic two-series negLnL scan overlay -- used for normalization-mode, no-singles, and dome-cut comparisons alike. seriesA is blue, seriesB red
def plot_neg_ln_l_scan_comparison(gridKev, negLnLCurveA, bestThresholdA, bestNegLnLA, labelA,
                                   negLnLCurveB, bestThresholdB, bestNegLnLB, labelB, savepath):
    plt.figure(figsize=(8, 6))
    plt.plot(gridKev, negLnLCurveA, 'o', markersize=3, color="steelblue", label=labelA)
    plt.axvline(bestThresholdA, color='steelblue', linestyle='--',
                label=f"{labelA} best fit: {bestThresholdA:0.2f} keV (" + r"$-2\Delta\ln L$" + f"={bestNegLnLA:0.2f})")
    plt.plot(gridKev, negLnLCurveB, 'o', markersize=3, color="red", label=labelB)
    plt.axvline(bestThresholdB, color='red', linestyle='--',
                label=f"{labelB} best fit: {bestThresholdB:0.2f} keV (" + r"$-2\Delta\ln L$" + f"={bestNegLnLB:0.2f})")
    plt.xlabel("Threshold [keV]", fontsize=16)
    plt.ylabel(r"$-2\Delta\ln L$", fontsize=16)
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.savefig(savepath)
    plt.close()

def plot_fit_to_seitz_ratio(groups, savepath, bestFitKey="bestFit"):
    seitzVals = [g["seitz"] for g in groups]
    fitVals = [g[bestFitKey]["threshold"] for g in groups]
    fitErrLow = [g[bestFitKey]["thresholdErrLow"] for g in groups]
    fitErrHigh = [g[bestFitKey]["thresholdErrHigh"] for g in groups]
    plt.figure(figsize=(8, 6))
    plt.errorbar(seitzVals, fitVals, yerr=[fitErrLow, fitErrHigh],
                 fmt='o', color="steelblue", ecolor="steelblue", capsize=4)
    plt.axline((0, 0), slope=1, color='gray', linestyle='--', label="Fit = Seitz")
    plt.xlabel("Seitz Threshold [keV]", fontsize=16)
    plt.ylabel("Best Fit Threshold [keV]", fontsize=16)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig(savepath)
    plt.close()

# generic two-series fit-vs-seitz overlay, same idea as plot_neg_ln_l_scan_comparison() above
def plot_fit_to_seitz_ratio_comparison(seitzVals, fitValsA, errLowA, errHighA, labelA,
                                        fitValsB, errLowB, errHighB, labelB, savepath, title=None):
    plt.figure(figsize=(8, 6))
    plt.errorbar(seitzVals, fitValsA, yerr=[errLowA, errHighA],
                 fmt='o', color="steelblue", ecolor="steelblue", capsize=4, label=labelA)
    plt.errorbar(seitzVals, fitValsB, yerr=[errLowB, errHighB],
                 fmt='o', color="red", ecolor="red", capsize=4, label=labelB)
    plt.axline((0, 0), slope=1, color='gray', linestyle='--', label="Fit = Seitz")
    plt.xlabel("Seitz Threshold [keV]", fontsize=16)
    plt.ylabel("Best Fit Threshold [keV]", fontsize=16)
    if title:
        plt.title(title, fontsize=14)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig(savepath)
    plt.close()

# pulls (threshold, thresholdErrLow, thresholdErrHigh) series out of a bestFit-shaped dict, for feeding plot_fit_to_seitz_ratio_comparison() above
def fit_threshold_series(groups, bestFitKey, thresholdField="threshold",
                          errLowField="thresholdErrLow", errHighField="thresholdErrHigh"):
    return (
        [g[bestFitKey][thresholdField] for g in groups],
        [g[bestFitKey][errLowField] for g in groups],
        [g[bestFitKey][errHighField] for g in groups],
    )

"""
STANDALONE COMBINED MULTIPLICITY PLOTS -- -2*deltaLnL of sim (at each group's own Seitz threshold,
no fitting) vs data, per group, for the standalone (single-series) plots skipped further up
"""
# The three pipelines differ only in which groups they hold, whether the sim side had the dome cut
# applied, and how they are labelled -- the statistics and both figures are identical work, so this is one
# loop rather than three copies of it. Stats are kept per pipeline because the final paper figure reuses
# the dome-cut ones. globalNormFactor is passed through but goes unused: every call here is at a group's
# own Seitz Q, where neg_ln_l_calc()/chi_squared_diagnostic() read g["seitzRate"] directly instead.
COMBINED_PIPELINES = [
    # name (drives the filename), groups, sim-side cut, corner label, chi2Debug.txt banner
    ("WithDomeCut",    groupsWithDomeCut,    ["dome"], "w/ dome cut",   "with dome cut"),
    ("WithoutDomeCut", groupsWithoutDomeCut, [],       "w/o dome cut",  "without dome cut"),
    ("NoCutAtAll",     groupsNoCutAtAll,     [],       "no cut at all", "no cut at all"),
]

combinedStats = {}
for pipelineName, pipelineGroups, simCut, cornerLabel, debugLabel in COMBINED_PIPELINES:
    globalNormFactor = global_normalization_factor(pipelineGroups)
    negLnLByGroup = [neg_ln_l_calc(g, g["seitz"], globalNormFactor, simCut) for g in pipelineGroups]
    chi2_debug(f"[chi2 debug] -- {debugLabel} --")
    chi2ByGroup = [chi_squared_diagnostic(g, g["seitz"], globalNormFactor, simCut, debug=True)
                   for g in pipelineGroups]
    combinedStats[pipelineName] = (negLnLByGroup, chi2ByGroup)

    for fileSuffix, showPerGroupStats in [("", False), ("PerGroupStats", True)]:
        plot_combined_multiplicity_single(
            pipelineGroups,
            savepath=output_path("comparison", f"combinedMultiplicity{pipelineName}{fileSuffix}.png"),
            negLnLByGroup=negLnLByGroup, chi2ByGroup=chi2ByGroup,
            sourceLabel=r'$^{252}$Cf Source:' + f'\n{cornerLabel}', showPerGroupStats=showPerGroupStats,
        )

# final paper version -- dome cut, no "w/ dome cut" line in the corner label, and the two highest
# threshold groups (1.91 and 3.10 keV) dropped; groups are sorted by seitz, so that's the tail
finalNegLnL, finalChi2 = combinedStats["WithDomeCut"]
plot_combined_multiplicity_single(
    groupsWithDomeCut[FINAL_PLOT_SLICE], savepath=output_path("comparison", "combinedMultiplicityFinal.pdf"),
    negLnLByGroup=finalNegLnL[FINAL_PLOT_SLICE], chi2ByGroup=finalChi2[FINAL_PLOT_SLICE],
    sourceLabel=r'$^{252}$Cf Source', groupsPerRow=3,
)

"""
SINGLE THRESHOLD RATES WITH AND WITHOUT THEORETICAL THRESHOLD
"""
binLabelsFull = ["1", "2", "3", "4", "5+"]

def plot_linhist(binLabels, binCounts, binCountErrorLow, binCountErrorHigh, backBins, backErrorLow, backErrorHigh,
                  backSubBins, backSubErrorLow, backSubErrorHigh, zeroKevRate, seitzRate, seitz,
                  savepath, bestFitRate=None, bestThreshold=None, bestNegLnL=None):
    plt.figure(figsize=(10, 10))
    x = np.arange(len(binLabels))
    edges = np.concatenate(([x[0] - 0.5], (x[:-1] + x[1:]) / 2, [x[-1] + 0.5]))

    plt.errorbar(x, binCounts, yerr=[binCountErrorLow, binCountErrorHigh], fmt='o', color="red", ecolor="red", label="Source Rate")
    plt.errorbar(x, backBins, yerr=[backErrorLow, backErrorHigh], fmt='o', color="blue", label="Background Rate")
    plt.errorbar(x, backSubBins, yerr=[backSubErrorLow, backSubErrorHigh], fmt='o', color="purple", label="Background Subtracted Rate")

    plt.stairs(zeroKevRate, edges, color="orange", linewidth=4, label="0keV Threshold")
    plt.stairs(seitzRate, edges, color="green", linewidth=6, label=f"Seitz Threshold\n({seitz:0.2f} keV )")
    if bestFitRate is not None:
        plt.stairs(bestFitRate, edges, color="magenta", linewidth=4, linestyle="--",
                   label=f"Best Fit Threshold\n({bestThreshold:0.2f} keV, " + r"$-2\Delta\ln L$" + f"={bestNegLnL:0.2f})")

    plt.xticks(x, binLabels, fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlabel("Bubble Multiplicity", fontsize=20)
    plt.ylabel("Rate [count/min]", fontsize=20)
    plt.legend(fontsize=20)
    plt.tight_layout()
    plt.savefig(savepath)
    plt.close()

"""
FULL PER-SIDE PIPELINE -- Z distributions, profile-likelihood fits, linhists, avg group. One function per stage; run_pipeline() at the bottom calls them in order, once per side.
"""

def plot_group_z_distributions(groups, outputDir):
    for g in groups[:numLowestSeitzZDist]:
        sourceZ = [z for pT in g["members"] for z in singleBubbleZByGroup.get(pT, [])]
        backgroundZ = [z for pT in g["members"] for z in backgroundZByGroup.get(pT, [])]
        plot_z_distribution(
            sourceZ, backgroundZ, g["seitz"],
            savepath=output_path(outputDir, f"zDistributions/zdist{g['label']}.png"),
            sourceLiveTime=g["liveTime"], backgroundLiveTime=g["backgroundTimeSec"] / 60,
        )


# fills every group's g["bestFit"] (normal, all-bins 1/2/3+ fit); if compareNormalizationModes is on, also fills g["bestFit"]["alt*"] for the side-by-side comparison plots only
def fit_groups(groups, globalNormFactor, simExcludedRegions):
    for g in groups:
        negLnLCurve = [neg_ln_l_calc(g, threshold, globalNormFactor, simExcludedRegions)
                       for threshold in thresholdRange]
        g["bestFit"] = compute_best_fit(g, negLnLCurve, simExcludedRegions)

        if compareNormalizationModes:
            altNegLnLCurve = [neg_ln_l_calc(g, threshold, globalNormFactor, simExcludedRegions,
                                             useGlobalNorm=not useGlobalNormalization)
                               for threshold in thresholdRange]
            altFit = compute_best_fit(g, altNegLnLCurve, simExcludedRegions)
            g["bestFit"].update({
                "altNegLnLCurve": altFit["negLnLCurve"],
                "altThreshold": altFit["threshold"],
                "altThresholdErrLow": altFit["thresholdErrLow"],
                "altThresholdErrHigh": altFit["thresholdErrHigh"],
                "altNegLnL": altFit["negLnL"],
            })


def plot_group_neg_ln_l_scans(groups, outputDir):
    for g in groups:
        plot_neg_ln_l_scan(
            g["bestFit"]["thresholdRange"], g["bestFit"]["negLnLCurve"],
            g["bestFit"]["threshold"], g["bestFit"]["negLnL"],
            savepath=output_path(outputDir, f"negLnLScans/negLnLScan{g['label']}.png"),
        )
        if compareNormalizationModes:
            plot_neg_ln_l_scan_comparison(
                g["bestFit"]["thresholdRange"], g["bestFit"]["negLnLCurve"],
                g["bestFit"]["threshold"], g["bestFit"]["negLnL"], normalization_mode_label(useGlobalNormalization),
                g["bestFit"]["altNegLnLCurve"], g["bestFit"]["altThreshold"], g["bestFit"]["altNegLnL"],
                normalization_mode_label(not useGlobalNormalization),
                savepath=output_path(outputDir, f"negLnLScans/negLnLScanCompare{g['label']}.png"),
            )


# separate fitting path, drops the multiplicity==1 bin -- stored in g["bestFitNoSingles"], never touches g["bestFit"] from fit_groups() above
def fit_groups_no_singles(groups, simExcludedRegions):
    for g in groups:
        negLnLCurveNoSingles = [neg_ln_l_calc_no_singles(g, threshold, simExcludedRegions)
                               for threshold in thresholdRange]
        g["bestFitNoSingles"] = compute_best_fit(g, negLnLCurveNoSingles, simExcludedRegions)


def plot_group_no_singles_comparisons(groups, outputDir):
    for g in groups:
        plot_neg_ln_l_scan(
            g["bestFitNoSingles"]["thresholdRange"], g["bestFitNoSingles"]["negLnLCurve"],
            g["bestFitNoSingles"]["threshold"], g["bestFitNoSingles"]["negLnL"],
            savepath=output_path(outputDir, f"negLnLScansNoSingles/negLnLScanNoSingles{g['label']}.png"),
        )
        plot_neg_ln_l_scan_comparison(
            g["bestFit"]["thresholdRange"], g["bestFit"]["negLnLCurve"],
            g["bestFit"]["threshold"], g["bestFit"]["negLnL"], "All bins (1, 2, 3+)",
            g["bestFitNoSingles"]["negLnLCurve"], g["bestFitNoSingles"]["threshold"], g["bestFitNoSingles"]["negLnL"],
            "No singles (2, 3+ only)",
            savepath=output_path(outputDir, f"negLnLScansNoSingles/negLnLScanCompare{g['label']}.png"),
        )

    plot_fit_to_seitz_ratio(
        groups, savepath=output_path(outputDir, "fitToSeitzRatioNoSingles.png"), bestFitKey="bestFitNoSingles"
    )
    plot_fit_to_seitz_ratio_comparison(
        [g["seitz"] for g in groups],
        *fit_threshold_series(groups, "bestFit"), "All bins (1, 2, 3+)",
        *fit_threshold_series(groups, "bestFitNoSingles"), "No singles (2, 3+ only)",
        savepath=output_path(outputDir, "fitToSeitzRatioCompareSingles.png"),
    )


# 0keV reference-line rate (no energy cut at all), rebinned and full versions -- scaled by the same
# normalizationFactor every group's own g["seitzRate"] uses (build_groups()'s single shared sim-to-data
# calibration), so 0keV, Seitz, and data sit on one consistent absolute scale across every panel instead
# of each panel independently re-normalizing 0keV to its own (noisy) observed total
def compute_zero_kev_reference(simExcludedRegions, normalizationFactor):
    zeroKevCountsRaw, _ = get_multiplicity_counts(0.0, excludedRegionsOverride=simExcludedRegions)
    zeroKevRateRebinned = [normalizationFactor * c for c in rebin(zeroKevCountsRaw)]
    zeroKevRateFull = [normalizationFactor * c for c in zeroKevCountsRaw]
    return zeroKevRateRebinned, zeroKevRateFull


# builds the positional args plot_linhist() expects for one group, respecting useRebinnedThresholdPlots --
# zeroKevRate is already globally scaled (see compute_zero_kev_reference above); the full/non-rebinned
# Seitz curve is put on that same global scale here too, matching how g["seitzRate"] already works in
# the rebinned branch -- shared by plot_group_linhists() and fit_seitz_multiplier()'s linHistsA plots
def build_linhist_args(g, zeroKevRateRebinned, zeroKevRateFull, normalizationFactor):
    if useRebinnedThresholdPlots:
        return (
            rebinnedBinLabels, g["binCounts"], g["binCountErrorLow"], g["binCountErrorHigh"],
            g["backBins"], g["backErrorLow"], g["backErrorHigh"],
            g["backSub"], g["errLow"], g["errHigh"],
            zeroKevRateRebinned, g["seitzRate"], g["seitz"],
        )
    seitzRateFull = [normalizationFactor * c for c in g["seitzCountsFull"]]
    return (
        binLabelsFull, g["binCountsFull"], g["binCountErrorLowFull"], g["binCountErrorHighFull"],
        g["backBinsFull"], g["backErrorLowFull"], g["backErrorHighFull"],
        g["backSubFull"], g["backSubErrorLowFull"], g["backSubErrorHighFull"],
        zeroKevRateFull, seitzRateFull, g["seitz"],
    )


def plot_group_linhists(groups, zeroKevRateRebinned, zeroKevRateFull, normalizationFactor, outputDir):
    for g in groups:
        linhistArgs = build_linhist_args(g, zeroKevRateRebinned, zeroKevRateFull, normalizationFactor)

        # one version without a best-fit curve, one with the normal (all-bins) fit, one with the no-singles (2, 3+ only) fit
        plot_linhist(*linhistArgs, savepath=output_path(outputDir, f"linHists/linhist{g['label']}.png"))
        if not computeBestFit:
            # no fitted threshold to overlay -- the plain linhist above is the whole story
            continue
        bestFitRate = g["bestFit"]["rate"] if useRebinnedThresholdPlots else g["bestFit"]["rateFull"]
        plot_linhist(
            *linhistArgs, savepath=output_path(outputDir, f"linHists/linhistFit{g['label']}.png"),
            bestFitRate=bestFitRate, bestThreshold=g["bestFit"]["threshold"], bestNegLnL=g["bestFit"]["negLnL"],
        )
        if computeNoSinglesFit:
            bestFitRateNoSingles = g["bestFitNoSingles"]["rate"] if useRebinnedThresholdPlots else g["bestFitNoSingles"]["rateFull"]
            plot_linhist(
                *linhistArgs, savepath=output_path(outputDir, f"linHists/linhistFitNoSingles{g['label']}.png"),
                bestFitRate=bestFitRateNoSingles, bestThreshold=g["bestFitNoSingles"]["threshold"],
                bestNegLnL=g["bestFitNoSingles"]["negLnL"],
            )


# the "average across every (p,T) group" version of the pipeline above: one combined rate group, one normal fit, one no-singles fit, three linhist variants
def run_avg_group_pipeline(groups, globalNormFactor, simExcludedRegions, positionDomeFlags,
                            backgroundPositionDomeFlags, zeroKevRateRebinned, zeroKevRateFull, outputDir):
    is_region_excluded = make_is_region_excluded(positionDomeFlags)
    background_is_region_excluded = make_is_region_excluded(backgroundPositionDomeFlags)
    # the avg group pools every (p, T)'s source events, so it pools their backgrounds too
    avgRateGroup = compute_rate_group(inGroup=lambda i: True, isVetoed=is_region_excluded,
                                      backgroundPTs=pToUse, backgroundIsVetoed=background_is_region_excluded)
    avgGroup = {
        **avgRateGroup,
        "errLow": avgRateGroup["backSubErrorLow"],
        "errHigh": avgRateGroup["backSubErrorHigh"],
    }

    avgSeitz = np.mean([g["seitz"] for g in groups])
    avgSeitzCountsRaw, _ = get_multiplicity_counts(avgSeitz, excludedRegionsOverride=simExcludedRegions)

    if computeBestFit:
        avgNegLnLCurve = [neg_ln_l_calc(avgGroup, threshold, globalNormFactor, simExcludedRegions)
                           for threshold in thresholdRange]
        avgBestFit = compute_best_fit(avgGroup, avgNegLnLCurve, simExcludedRegions)

        if compareNormalizationModes:
            avgAltNegLnLCurve = [neg_ln_l_calc(avgGroup, threshold, globalNormFactor, simExcludedRegions,
                                                useGlobalNorm=not useGlobalNormalization)
                                  for threshold in thresholdRange]
            avgAltFit = compute_best_fit(avgGroup, avgAltNegLnLCurve, simExcludedRegions)
            avgBestFit.update({
                "altNegLnLCurve": avgAltFit["negLnLCurve"],
                "altThreshold": avgAltFit["threshold"],
                "altNegLnL": avgAltFit["negLnL"],
            })

        if computeNoSinglesFit:
            # separate no-singles fit for the avg group, same as fit_groups_no_singles() above
            avgNegLnLCurveNoSingles = [neg_ln_l_calc_no_singles(avgGroup, threshold, simExcludedRegions)
                                      for threshold in thresholdRange]
            avgBestFitNoSingles = compute_best_fit(avgGroup, avgNegLnLCurveNoSingles, simExcludedRegions)

            plot_neg_ln_l_scan(
                avgBestFitNoSingles["thresholdRange"], avgBestFitNoSingles["negLnLCurve"],
                avgBestFitNoSingles["threshold"], avgBestFitNoSingles["negLnL"],
                savepath=output_path(outputDir, "negLnLScanAvgNoSingles.png"),
            )
            plot_neg_ln_l_scan_comparison(
                avgBestFit["thresholdRange"], avgBestFit["negLnLCurve"],
                avgBestFit["threshold"], avgBestFit["negLnL"], "All bins (1, 2, 3+)",
                avgBestFitNoSingles["negLnLCurve"], avgBestFitNoSingles["threshold"], avgBestFitNoSingles["negLnL"],
                "No singles (2, 3+ only)",
                savepath=output_path(outputDir, "negLnLScanAvgCompareSingles.png"),
            )

        plot_neg_ln_l_scan(
            avgBestFit["thresholdRange"], avgBestFit["negLnLCurve"],
            avgBestFit["threshold"], avgBestFit["negLnL"],
            savepath=output_path(outputDir, "negLnLScanAvg.png"),
        )
        if compareNormalizationModes:
            plot_neg_ln_l_scan_comparison(
                avgBestFit["thresholdRange"], avgBestFit["negLnLCurve"],
                avgBestFit["threshold"], avgBestFit["negLnL"], normalization_mode_label(useGlobalNormalization),
                avgBestFit["altNegLnLCurve"], avgBestFit["altThreshold"], avgBestFit["altNegLnL"],
                normalization_mode_label(not useGlobalNormalization),
                savepath=output_path(outputDir, "negLnLScanCompareAvg.png"),
            )

    if useRebinnedThresholdPlots:
        avgSeitzCountsRebinned = rebin(avgSeitzCountsRaw)
        # avgSeitzRate stays on its own per-total scale -- normalizationFactor isn't reliable here since
        # avgSeitz usually falls outside the thresholds it was calibrated on. zeroKevRateRebinned is
        # already globally scaled (see compute_zero_kev_reference), so it's used as-is, unlike avgSeitzRate
        avgSeitzRate = seitz_count(counts_to_ratios(avgSeitzCountsRebinned), sum(avgGroup["backSubUnfloored"]))
        avgLinhistArgs = (
            rebinnedBinLabels, avgGroup["binCounts"], avgGroup["binCountErrorLow"], avgGroup["binCountErrorHigh"],
            avgGroup["backBins"], avgGroup["backErrorLow"], avgGroup["backErrorHigh"],
            avgGroup["backSub"], avgGroup["errLow"], avgGroup["errHigh"],
            zeroKevRateRebinned, avgSeitzRate, avgSeitz,
        )
    else:
        totalAvg = sum(avgGroup["backSubUnflooredFull"])
        avgSeitzRatiosFull = counts_to_ratios(avgSeitzCountsRaw)
        avgLinhistArgs = (
            binLabelsFull, avgGroup["binCountsFull"], avgGroup["binCountErrorLowFull"], avgGroup["binCountErrorHighFull"],
            avgGroup["backBinsFull"], avgGroup["backErrorLowFull"], avgGroup["backErrorHighFull"],
            avgGroup["backSubFull"], avgGroup["backSubErrorLowFull"], avgGroup["backSubErrorHighFull"],
            zeroKevRateFull, seitz_count(avgSeitzRatiosFull, totalAvg), avgSeitz,
        )

    # one version without a best-fit curve, one with the normal (all-bins) fit, one with the no-singles (2, 3+ only) fit
    plot_linhist(*avgLinhistArgs, savepath=output_path(outputDir, "avgseitz.png"))
    if computeBestFit:
        avgBestFitRate = avgBestFit["rate"] if useRebinnedThresholdPlots else avgBestFit["rateFull"]
        plot_linhist(
            *avgLinhistArgs, savepath=output_path(outputDir, "avgseitzFit.png"),
            bestFitRate=avgBestFitRate, bestThreshold=avgBestFit["threshold"], bestNegLnL=avgBestFit["negLnL"],
        )
        if computeNoSinglesFit:
            avgBestFitRateNoSingles = avgBestFitNoSingles["rate"] if useRebinnedThresholdPlots else avgBestFitNoSingles["rateFull"]
            plot_linhist(
                *avgLinhistArgs, savepath=output_path(outputDir, "avgseitzFitNoSingles.png"),
                bestFitRate=avgBestFitRateNoSingles, bestThreshold=avgBestFitNoSingles["threshold"],
                bestNegLnL=avgBestFitNoSingles["negLnL"],
            )


def run_pipeline(groups, normalizationFactor, simExcludedRegions, positionDomeFlags,
                 backgroundPositionDomeFlags, outputDir):
    plot_group_z_distributions(groups, outputDir)

    globalNormFactor = global_normalization_factor(groups)
    if computeBestFit:
        fit_groups(groups, globalNormFactor, simExcludedRegions)
        plot_group_neg_ln_l_scans(groups, outputDir)

        plot_fit_to_seitz_ratio(groups, savepath=output_path(outputDir, "fitToSeitzRatio.png"))
        if compareNormalizationModes:
            plot_fit_to_seitz_ratio_comparison(
                [g["seitz"] for g in groups],
                *fit_threshold_series(groups, "bestFit"), normalization_mode_label(useGlobalNormalization),
                *fit_threshold_series(groups, "bestFit", "altThreshold", "altThresholdErrLow", "altThresholdErrHigh"),
                normalization_mode_label(not useGlobalNormalization),
                savepath=output_path(outputDir, "fitToSeitzRatioCompare.png"),
            )

        if computeNoSinglesFit:
            fit_groups_no_singles(groups, simExcludedRegions)
            plot_group_no_singles_comparisons(groups, outputDir)

    zeroKevRateRebinned, zeroKevRateFull = compute_zero_kev_reference(simExcludedRegions, normalizationFactor)
    plot_group_linhists(groups, zeroKevRateRebinned, zeroKevRateFull, normalizationFactor, outputDir)

    run_avg_group_pipeline(groups, globalNormFactor, simExcludedRegions, positionDomeFlags,
                            backgroundPositionDomeFlags, zeroKevRateRebinned, zeroKevRateFull, outputDir)


# total negLnL across every group when every group is fit against threshold = A * g["seitz"], for one shared A instead of a separate best-fit threshold per group
def neg_ln_l_calc_seitz_multiplier(groups, A, globalNormFactor, simExcludedRegions):
    return sum(neg_ln_l_calc(g, A * g["seitz"], globalNormFactor, simExcludedRegions) for g in groups)

# chi2 analog of neg_ln_l_calc_seitz_multiplier() above, diagnostic-only (see chi_squared_diagnostic())
def chi_squared_diagnostic_seitz_multiplier(groups, A, globalNormFactor, simExcludedRegions):
    return sum(chi_squared_diagnostic(g, A * g["seitz"], globalNormFactor, simExcludedRegions) for g in groups)

# overlays negLnL and chi2 on the same axis/scale -- under Wilks' theorem they should track each other
# reasonably closely, so this is a sanity check that the profile-likelihood swap didn't produce something
# qualitatively different from the error-based chi2 this pipeline used before
def plot_neg_ln_l_vs_chi2_diagnostic(gridA, negLnLCurve, chi2Curve, savepath, xlabel="Seitz threshold multiplier A"):
    plt.figure(figsize=(8, 6))
    plt.plot(gridA, negLnLCurve, 'o', markersize=3, color="steelblue", label=r"$-2\Delta\ln L$")
    plt.plot(gridA, chi2Curve, 'o', markersize=3, color="darkorange", label=r"$\chi^2$")
    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel("Test statistic", fontsize=16)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig(savepath)
    plt.close()


# single-parameter fit: scans A over seitzMultiplierRange, picks the A minimizing the combined negLnL across every (p,T) group at once, then reports it per group as threshold = A * g["seitz"] (stored in g["bestFitA"], alongside g["bestFit"]/g["bestFitNoSingles"]). Run once per side, after run_pipeline(). All plots land in outputDir/linHistsA/
def fit_seitz_multiplier(groups, simExcludedRegions, outputDir, normalizationFactor):
    globalNormFactor = global_normalization_factor(groups)
    negLnLCurve = [neg_ln_l_calc_seitz_multiplier(groups, A, globalNormFactor, simExcludedRegions)
                 for A in seitzMultiplierRange]
    bestIdx = int(np.argmin(negLnLCurve))
    bestA = seitzMultiplierRange[bestIdx]
    lowA, highA = neg_ln_l_confidence_interval(seitzMultiplierRange, negLnLCurve, bestIdx)
    bestNegLnL = negLnLCurve[bestIdx]
    AErrLow, AErrHigh = bestA - lowA, highA - bestA

    # the rebinned (1, 2, 3+ or 1, 2+) or 5 full (1, 2, 3, 4, 5+) multiplicity bins per (p, T) group
    binsPerGroup = len(rebinnedBinLabels) if useRebinnedThresholdPlots else 5
    nBins = len(groups) * binsPerGroup

    plot_neg_ln_l_scan(
        seitzMultiplierRange, negLnLCurve, bestA, bestNegLnL,
        savepath=output_path(outputDir, "negLnLScanSeitzMultiplier.png"),
        xlabel="Seitz threshold multiplier A",
        bestLabel=f"Best fit: A = {bestA:0.3f} (" + r"$-2\Delta\ln L$" + f"={bestNegLnL:0.2f})",
        sigmaRange=(lowA, highA),
        infoText=f"N bins = {nBins} ({len(groups)} (p, T) groups x {binsPerGroup} multiplicity bins)",
    )

    chi2Curve = [chi_squared_diagnostic_seitz_multiplier(groups, A, globalNormFactor, simExcludedRegions)
                 for A in seitzMultiplierRange]
    plot_neg_ln_l_vs_chi2_diagnostic(
        seitzMultiplierRange, negLnLCurve, chi2Curve,
        savepath=output_path(outputDir, "negLnLScanSeitzMultiplierChi2Diagnostic.png"),
    )

    # per-group threshold = A * seitz, with its own negLnL and A's uncertainty propagated multiplicatively
    for g in groups:
        threshold = bestA * g["seitz"]
        groupNegLnL = neg_ln_l_calc(g, threshold, globalNormFactor, simExcludedRegions)
        g["bestFitA"] = build_fit_result(
            g, threshold, groupNegLnL, simExcludedRegions,
            thresholdErrLow=AErrLow * g["seitz"], thresholdErrHigh=AErrHigh * g["seitz"],
        )

    # same trio of plots the normal fit gets (per-group negLnL scan, fit-vs-seitz scatter, linhist with the fit overlaid), all under linHistsA/
    plot_fit_to_seitz_ratio(groups, savepath=output_path(outputDir, "linHistsA/fitToSeitzRatioA.png"), bestFitKey="bestFitA")

    zeroKevRateRebinned, zeroKevRateFull = compute_zero_kev_reference(simExcludedRegions, normalizationFactor)
    for g in groups:
        plot_neg_ln_l_scan(
            g["bestFit"]["thresholdRange"], g["bestFit"]["negLnLCurve"],
            g["bestFitA"]["threshold"], g["bestFitA"]["negLnL"],
            savepath=output_path(outputDir, f"linHistsA/negLnLScanA{g['label']}.png"),
            bestLabel=f"A * Seitz fit: {g['bestFitA']['threshold']:0.2f} keV (" + r"$-2\Delta\ln L$" + f"={g['bestFitA']['negLnL']:0.2f})",
        )

        linhistArgs = build_linhist_args(g, zeroKevRateRebinned, zeroKevRateFull, normalizationFactor)
        bestFitRateA = g["bestFitA"]["rate"] if useRebinnedThresholdPlots else g["bestFitA"]["rateFull"]
        plot_linhist(
            *linhistArgs, savepath=output_path(outputDir, f"linHistsA/linhistFitA{g['label']}.png"),
            bestFitRate=bestFitRateA, bestThreshold=g["bestFitA"]["threshold"], bestNegLnL=g["bestFitA"]["negLnL"],
        )

    return {
        "aRange": seitzMultiplierRange,
        "negLnLCurve": negLnLCurve,
        "A": bestA,
        "AErrLow": AErrLow,
        "AErrHigh": AErrHigh,
        "negLnL": bestNegLnL,
    }


run_pipeline(groupsWithoutDomeCut, normalizationFactorWithoutDomeCut, [], positionDomeFlagsPre,
            backgroundPositionDomeFlagsPre, "withoutDomeCut")
run_pipeline(groupsWithDomeCut, normalizationFactorWithDomeCut, ["dome"], positionDomeFlagsPost,
            backgroundPositionDomeFlagsPost, "withDomeCut")

# everything below is fit output -- skipped wholesale when computeBestFit is off
if computeBestFit:
    # dome cut vs no dome cut, normal (all-bins, 1/2/3+) fit -- needs both run_pipeline() calls above to have already populated g["bestFit"]
    plot_fit_to_seitz_ratio_comparison(
        [g["seitz"] for g in groupsWithDomeCut],
        *fit_threshold_series(groupsWithDomeCut, "bestFit"), "Dome cut",
        *fit_threshold_series(groupsWithoutDomeCut, "bestFit"), "No dome cut",
        savepath=output_path("comparison", "fitToSeitzRatioDomeCutCompare.png"),
        title="All bins (1, 2, 3+) fit",
    )

    # same comparison, broken out per (p,T) group as negLnL-vs-threshold scans
    for gWithCut, gWithoutCut in zip(groupsWithDomeCut, groupsWithoutDomeCut):
        plot_neg_ln_l_scan_comparison(
            gWithCut["bestFit"]["thresholdRange"], gWithCut["bestFit"]["negLnLCurve"],
            gWithCut["bestFit"]["threshold"], gWithCut["bestFit"]["negLnL"], "Dome cut",
            gWithoutCut["bestFit"]["negLnLCurve"], gWithoutCut["bestFit"]["threshold"], gWithoutCut["bestFit"]["negLnL"],
            "No dome cut",
            savepath=output_path("comparison", f"negLnLScans/negLnLScanDomeCutCompare{gWithCut['label']}.png"),
        )

    # single-A "threshold = A * seitz" fit, run last: once for pre-cut (withoutDomeCut), once for post-cut (withDomeCut)
    fit_seitz_multiplier(groupsWithoutDomeCut, [], "withoutDomeCut", normalizationFactorWithoutDomeCut)
    fit_seitz_multiplier(groupsWithDomeCut, ["dome"], "withDomeCut", normalizationFactorWithDomeCut)
