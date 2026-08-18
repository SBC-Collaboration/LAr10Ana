## Variant of neutronExcluded.py: identical pipeline, except build_groups()'s normalizationFactor is
## livetime-weighted instead of a plain per-group mean (see the comment at that computation below).
## Quick A/B comparison copy, not meant to replace neutronExcluded.py -- writes to its own PLOTS_ROOT
## so it never overwrites the main pipeline's output.
## dome exclusion: same rate/threshold-fit pipeline as neutronPlotting.py, but the sim and real sides both drop dome-region hits, gated by excludedRegions below
import atexit
import glob
import os
import sys

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
RECON_DIR = "/exp/e961/data/SBC-25-recon/v0.4.2/"

DOME_Z_THRESHOLD_CM = -3

SIM_DOME_Z_THRESHOLD_MM = 591

# drip cut then dome cut
DOME_Z_THRESHOLD_CM_PRE = 6
DOME_Z_THRESHOLD_CM_POST = -3


PLOTS_ROOT = "plots"


def output_path(outputDir, filename):
    path = os.path.join(PLOTS_ROOT, outputDir, filename)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path


# every chi_squared_diagnostic(debug=True) line goes here instead of the console; opened on the first write, truncated once per run
chi2DebugFile = None

def chi2_debug(line):
    global chi2DebugFile
    if chi2DebugFile is None:
        chi2DebugFile = open(output_path("comparison", "chi2Debug.txt"), "w")
        atexit.register(chi2DebugFile.close)
    chi2DebugFile.write(line + "\n")


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
backgroundRunsHot = ["20260217_7","20260217_8","20260217_9","20260217_10","20260217_11","20260217_12","20260217_13","20260218_0","20260218_1","20260218_2","20260218_3","20260218_4","20260218_5","20260218_6","20260218_15","20260218_16","20260219_0","20260219_1","20260219_2","20260219_3","20260219_4","20260219_5","20260219_6","20260219_7","20260219_8","20260219_9","20260219_10","20260219_11","20260220_1","20260220_2","20260220_3","20260220_4",]
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
# background events without a fixed setpoint belong to no group and are dropped, same rule pToUse applies to the source side
def load_background():
    subdirs = {os.path.basename(p.rstrip(os.sep)) for p in glob.glob(os.path.join(RECON_DIR, '*/'))}
    byPT = {}
    expDataCache, evDataCache = {}, {}
    for run, ev, mult, _ in iter_matched_events(HANDSCAN_DIR, backgroundList):
        if run not in subdirs:
            continue
        if run not in expDataCache:
            expDataCache[run] = Streamer(os.path.join(RECON_DIR, run, 'exposure.sbc')).to_dict()
        expData = expDataCache[run]
        for i in range(len(expData["ev"])):
            if int(expData["ev"][i]) == int(ev) and float(expData['PT2121_livetime'][i]) > 1:
                pset = _event_pset_temp(run, ev, evDataCache)
                if pset is not None and pset[0] == pset[1]:
                    binCounts, liveTime = byPT.setdefault((pset[0], pset[2]), ([0] * 5, 0.0))
                    if mult != 0:
                        binCounts[min(mult, 5) - 1] += 1
                    byPT[(pset[0], pset[2])] = (binCounts, liveTime + float(expData['PT2121_livetime'][i]))
                break
    return byPT

# (mult, region) per neutron-run event, its live time, and its (pset_lo, pset_hi, temp)
def load_neutron_events():
    bubbleCount, sourceTimes, psetsTemps, runEvs = [], [], [], []
    expDataCache, evDataCache = {}, {}
    for run, ev, mult, region in iter_matched_events(HANDSCAN_DIR, neutronRuns):
        if run not in expDataCache:
            expDataCache[run] = Streamer(f'{RECON_DIR}{run}/exposure.sbc').to_dict()
        expData = expDataCache[run]
        liveTime = None
        for i in range(len(expData["ev"])):
            if int(expData["ev"][i]) == int(ev):
                liveTime = float(expData['PT2121_livetime'][i])
                break
        pset = _event_pset_temp(run, ev, evDataCache)
        # all four lists are indexed in parallel downstream, so an event the exposure or event file doesn't
        # have goes in none of them -- appending to only some would shift every later event's livetime/(p, T)
        if liveTime is None or pset is None:
            continue
        bubbleCount.append((mult, region))
        runEvs.append((run, ev))
        sourceTimes.append(liveTime)
        psetsTemps.append(pset)
    return bubbleCount, sourceTimes, psetsTemps, runEvs


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

# counts per multiplicity bin (1..5+) and total live time, for events where keep(i, region) is True (for now always but later we could exclude wall/bellows/etc.)
def bin_multiplicities(bubbleCount, sourceTimes, keep):
    binCounts = [0] * 5
    sourceTime = 0.0
    for i, (mult, region) in enumerate(bubbleCount):
        if not keep(i, region):
            continue
        sourceTime += sourceTimes[i]
        if mult != 0:
            binCounts[min(mult, 5) - 1] += 1
    return binCounts, sourceTime


# background-subtracted counts and their (asymmetric) errors, background scaled to sourceTime
def background_subtract(binCounts, sourceTime, backgroundBinCounts, backgroundTime):
    scale = sourceTime / backgroundTime if backgroundTime > 0 else 0.0
    backBins = [b * scale for b in backgroundBinCounts]
    backErrorLow = [np.sqrt(b) * scale if b >= 1 else 0.0 for b in backgroundBinCounts]
    backErrorHigh = [np.sqrt(b) * scale if b >= 1 else -np.log(1 - 0.68) * scale for b in backgroundBinCounts]
    binCountError = [np.sqrt(c) for c in binCounts]
    backSubBins = [c - b for c, b in zip(binCounts, backBins)]
    # background high -> subtracted rate pulled down; background low -> subtracted rate pulled up
    backSubErrorLow = [np.sqrt(countErr**2 + backErr**2) for countErr, backErr in zip(binCountError, backErrorHigh)]
    backSubErrorHigh = [np.sqrt(countErr**2 + backErr**2) for countErr, backErr in zip(binCountError, backErrorLow)]
    return backBins, backErrorLow, backErrorHigh, backSubBins, backSubErrorLow, backSubErrorHigh, binCountError


# PLOTTING ONLY -- a negative rate is unphysical, so the drawn point sits at 0 with no lower error and an upper error
# trimmed to the 68% upper edge the fluctuation actually allows. Deliberately NOT fed to the fits: flooring an observed
# value pulls it toward the prediction, and the trimmed error becomes the chi2 denominator, so a deficit bin either
# blows up or (once its upper edge reaches 0) drops out of the sum entirely. Statistics use the unfloored values below
def floor_at_zero(values, errLow, errHigh):
    flooredErrHigh = [max(v + e, 0.0) if v < 0 else e for v, e in zip(values, errHigh)]
    flooredErrLow = [min(e, v) if v > 0 else 0.0 for v, e in zip(values, errLow)]
    return [max(v, 0.0) for v in values], flooredErrLow, flooredErrHigh


# collapse the 5 multiplicity classes (1,2,3,4,5+) into 3: [1, 2, 3+]
def rebin(values):
    return [values[0], values[1], sum(values[2:])]

# values are independent counts/rates, so combine their errors in quadrature
def rebin_errors(errors):
    return [errors[0], errors[1], np.sqrt(sum(e**2 for e in errors[2:]))]


def counts_to_ratios(counts):
    total = sum(counts)
    return [c / total for c in counts]

# scale a multiplicity ratio shape to match a target total
def seitz_count(ratios, total):
    scale = total / sum(ratios)
    return [scale * r for r in ratios]

# load in data (excludedRegions-independent -- region filtering happens per-group below)
backgroundByPT = load_background()

# background counts/live time for one (p, T), or summed over several for the all-groups-averaged pipeline
def background_for(pTs):
    binCounts = [0] * 5
    liveTime = 0.0
    for pt in pTs:
        counts, t = backgroundByPT.get(pt, ([0] * 5, 0.0))
        binCounts = [a + b for a, b in zip(binCounts, counts)]
        liveTime += t
    return binCounts, liveTime

bubbleCount, sourceTimes, psetsTemps, neutronRunEvs = load_neutron_events()
singleBubblePositions = load_single_bubble_positions(neutronRuns, label="source")

# reco-Z-based dome flags for the withoutDomeCut/withDomeCut pipeline: "pre" = always-applied loose threshold, "post" = tighter threshold on top
positionDomeFlagsPre = compute_position_dome_flags(bubbleCount, neutronRunEvs, DOME_Z_THRESHOLD_CM_PRE)
positionDomeFlagsPost = compute_position_dome_flags(bubbleCount, neutronRunEvs, DOME_Z_THRESHOLD_CM_POST)

# (p, T) pairs with a fixed pressure setpoint, deduped
pToUse = sorted({(float(lo), float(t)) for lo, hi, t in psetsTemps if float(lo) == float(hi)})


# bins events matching `keep`, background-subtracts, and converts to a rate in counts/minute -- shared core of both build_groups()'s per-(p,T) loop and the avg-group section
def compute_rate_group(keep, backgroundPTs):
    binCounts, liveTimeSec = bin_multiplicities(bubbleCount, sourceTimes, keep)
    backgroundBinCounts, backgroundTime = background_for(backgroundPTs)
    backBins, backErrorLow, backErrorHigh, backSubBins, backSubErrorLow, backSubErrorHigh, binCountError = \
        background_subtract(binCounts, liveTimeSec, backgroundBinCounts, backgroundTime)

    flooredBins, flooredErrorLow, flooredErrorHigh = floor_at_zero(backSubBins, backSubErrorLow, backSubErrorHigh)

    liveTimeMin = liveTimeSec / 60
    toRate = lambda values: [v / liveTimeMin for v in values]
    return {
        "liveTime": liveTimeMin,
        "liveTimeSec": liveTimeSec,
        # raw (un-scaled, un-rated) inputs to the subtraction, for the chi2 debug dump
        "binCountsRaw": binCounts,
        "backCountsRaw": list(backgroundBinCounts),
        "backgroundTimeSec": backgroundTime,
        "backgroundScale": liveTimeSec / backgroundTime if backgroundTime > 0 else 0.0,
        "binCounts": toRate(binCounts),
        "binCountError": toRate(binCountError),
        "backBins": toRate(backBins),
        "backErrorLow": toRate(backErrorLow),
        "backErrorHigh": toRate(backErrorHigh),
        # plotted (floored at 0) ...
        "backSubFull": toRate(flooredBins),
        "backSubErrorLowFull": toRate(flooredErrorLow),
        "backSubErrorHighFull": toRate(flooredErrorHigh),
        # ... and the same rates as actually subtracted, which the fits and chi2 use
        "backSubUnflooredFull": toRate(backSubBins),
        "errLowUnflooredFull": toRate(backSubErrorLow),
        "errHighUnflooredFull": toRate(backSubErrorHigh),
    }


# builds the full per-(p,T) "groups" list for one (simExcludedRegions, positionDomeFlags) combination -- lets the pre/post pipeline independently control the sim-side cut and real-side position cut
def build_groups(simExcludedRegions, positionDomeFlags):
    groups = []
    is_region_excluded = make_is_region_excluded(positionDomeFlags)

    for p, T in pToUse:
        rateGroup = compute_rate_group(
            keep=lambda i, region, p=p, T=T: not is_region_excluded(i, region)
                                              and psetsTemps[i][0] == p and psetsTemps[i][2] == T,
            backgroundPTs=[(p, T)],
        )

        # seitz threshold for this (P,T) pair, fed straight into the dome-excluded Cf sim counts
        seitz = sm.SeitzModel(p * 14.5038, -273.15 + T, 'argon').Q
        seitzCounts, _ = get_multiplicity_counts(seitz, excludedRegionsOverride=simExcludedRegions)
        groups.append({
            "p": p,
            "T": T,
            "seitz": seitz,
            "seitzCountsFull": seitzCounts,
            "backSub": rebin(rateGroup["backSubFull"]),
            "errLow": rebin_errors(rateGroup["backSubErrorLowFull"]),
            "errHigh": rebin_errors(rateGroup["backSubErrorHighFull"]),
            "backSubUnfloored": rebin(rateGroup["backSubUnflooredFull"]),
            "errLowUnfloored": rebin_errors(rateGroup["errLowUnflooredFull"]),
            "errHighUnfloored": rebin_errors(rateGroup["errHighUnflooredFull"]),
            "seitzRate": rebin(seitzCounts),
            "seitzCounts": rebin(seitzCounts),
            "bestFit": {},
            **rateGroup,
        })

    # sort from lowest to highest seitz threshold
    groups.sort(key=lambda g: g["seitz"])

    # scale sim-predicted counts to match observed data rates, weighted by each group's own livetime
    # instead of a plain mean -- N = sum_i(T_i * ratio_i) / sum_i(T_i), so groups with more exposure
    # (and thus a more reliable data/sim ratio) pull the shared factor harder than short, noisy ones
    observedRatesByGroup = [g["backSubUnfloored"] for g in groups]
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
groupsWithoutDomeCut, normalizationFactorWithoutDomeCut = build_groups([], positionDomeFlagsPre)
groupsWithDomeCut, normalizationFactorWithDomeCut = build_groups(["dome"], positionDomeFlagsPost)

# fully uncut baseline: sim off, real side has neither the PRE (always-applied) nor the POST dome
# threshold -- combined-plot-only, doesn't get the full run_pipeline() treatment (no linhists/scans/Z-dists)
positionDomeFlagsNone = [False] * len(bubbleCount)
groupsNoCutAtAll, normalizationFactorNoCutAtAll = build_groups([], positionDomeFlagsNone)

# dumps every group's background-subtracted data rate and Seitz/sim-predicted rate, pre-cut vs post-cut, to a plain text file
def write_rates_txt(groupsPre, groupsPost, savepath):
    binLabels = ["1", "2", "3+"]
    with open(savepath, "w") as f:
        f.write(f"Rates pre-cut (uncut baseline: real side Z > {DOME_Z_THRESHOLD_CM_PRE}cm, sim off) "
                f"vs post-cut (real side Z > {DOME_Z_THRESHOLD_CM_POST}cm, sim on), count/min\n")
        f.write("data = background-subtracted real rate, sim = Seitz/sim-predicted rate\n")
        f.write("=" * 78 + "\n\n")
        for gPre, gPost in zip(groupsPre, groupsPost):
            f.write(f"Seitz = {gPre['seitz']:0.2f} keV  (P = {gPre['p']:0.2f} bar, T = {gPre['T']:0.1f} K)\n")
            f.write(f"  Pre-cut  (Z > {DOME_Z_THRESHOLD_CM_PRE}cm, sim off):  livetime = {gPre['liveTime']:0.2f} min\n")
            for label, rate, errLow, errHigh, simRate in zip(
                    binLabels, gPre["backSub"], gPre["errLow"], gPre["errHigh"], gPre["seitzRate"]):
                f.write(f"    mult {label:<3}: data = {rate:0.4f} (+{errHigh:0.4f}/-{errLow:0.4f})   sim = {simRate:0.4f}\n")
            f.write(f"  Post-cut (Z > {DOME_Z_THRESHOLD_CM_POST}cm, sim on): livetime = {gPost['liveTime']:0.2f} min\n")
            for label, rate, errLow, errHigh, simRate in zip(
                    binLabels, gPost["backSub"], gPost["errLow"], gPost["errHigh"], gPost["seitzRate"]):
                f.write(f"    mult {label:<3}: data = {rate:0.4f} (+{errHigh:0.4f}/-{errLow:0.4f})   sim = {simRate:0.4f}\n")
            f.write("\n")

write_rates_txt(groupsWithoutDomeCut, groupsWithDomeCut, savepath=output_path("comparison", "ratesPrePostCut.txt"))

## plot making
"""
COMBINED PAPER PLOT -- side by side, dome cut vs no dome cut
"""
# groupsWithCut/groupsWithoutCut must be the same length, in the same (p,T) order (guaranteed since both come from build_groups() over the same pToUse, sorted by seitz)
def plot_combined_multiplicity_comparison(groupsWithCut, groupsWithoutCut, savepath, groupsPerRow=3):
    binLabels = ["1", "2", "3+"]
    numBins = 3
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
    binLabels = ["1", "2", "3+"]
    numBins = 3
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
                perGroupText = (r'$-2\Delta\ln\mathcal{L}$' + f'={groupNegLnL:0.1f}\n'
                                 + r'$\chi^2$' + f'={groupChi2:0.1f}')
                ax.text(center, 0.86, perGroupText, transform=trans, ha='center', va='top', fontsize=6)

            if gi < len(rowGroups) - 1:
                ax.axvline(pos + numBins + gap / 2 - 0.5, linestyle='--', linewidth=0.7,
                           color='gray', zorder=0)

            pos += numBins + gap

        ax.set_ylim(0, globalMax * 1.2)
        ax.set_xlim(-1, pos - gap)
        ax.set_xticks([])
        ax.tick_params(axis='y', labelsize=12)

    totalNegLnL = sum(negLnLByGroup)
    totalChi2 = sum(chi2ByGroup)
    totalBins = len(groups) * numBins
    statText = r'$-2\Delta\ln\mathcal{L}$/n.d.o.f.: ' + f'{totalNegLnL:0.1f}/({totalBins} - 1)'
    chi2Text = r'$\chi^2$/n.d.o.f.: ' + f'{totalChi2:0.1f}/({totalBins} - 1)'

    axes[0].set_title(r'$Q_{seitz}$ [keV]', loc='left', fontsize=16, pad=2)
    axes[0].set_title(sourceLabel, loc='right', fontsize=12, pad=2)
    axes[nRows // 2].set_ylabel("Rate [count/min]", fontsize=16)
    axes[-1].set_xlabel("Bubble Multiplicity", fontsize=12, labelpad=28)

    if showPerGroupStats:
        # per-group stat text already occupies the top of every panel (including the last row's), so the
        # combined total goes below the whole figure instead of inside axes[-1] -- putting it there too
        # overlapped the last panel's own per-group numbers
        fig.tight_layout()
        fig.subplots_adjust(hspace=0.15, bottom=0.16)
        fig.text(0.5, 0.02, statText + '     ' + chi2Text, ha='center', va='bottom', fontsize=9)
    else:
        axes[-1].text(0.99, 0.80, statText + '\n' + chi2Text, transform=axes[-1].transAxes,
                      ha='right', va='top', fontsize=10)
        fig.tight_layout()
        fig.subplots_adjust(hspace=0.15)

    fig.savefig(savepath)
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
        binLabels = ["1", "2", "3+"]
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
        pT = f'(p={dataGroup["p"]:0.2f}, T={dataGroup["T"]:0.1f})' if "p" in dataGroup else "(avg group)"
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

# same idea as neg_ln_l_calc, but drops the multiplicity==1 bin and fits only 2 and 3+ -- always rebinned, results go in g["bestFitNoSingles"], never g["bestFit"]
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
globalNormFactorWithDomeCut = global_normalization_factor(groupsWithDomeCut)
negLnLByGroupWithDomeCut = [
    neg_ln_l_calc(g, g["seitz"], globalNormFactorWithDomeCut, ["dome"]) for g in groupsWithDomeCut
]
chi2_debug("[chi2 debug] -- with dome cut --")
chi2ByGroupWithDomeCut = [
    chi_squared_diagnostic(g, g["seitz"], globalNormFactorWithDomeCut, ["dome"], debug=True)
    for g in groupsWithDomeCut
]
plot_combined_multiplicity_single(
    groupsWithDomeCut, savepath=output_path("comparison", "combinedMultiplicityWithDomeCut.png"),
    negLnLByGroup=negLnLByGroupWithDomeCut, chi2ByGroup=chi2ByGroupWithDomeCut,
    sourceLabel=r'$^{252}$Cf Source:' + '\nw/ dome cut',
)
plot_combined_multiplicity_single(
    groupsWithDomeCut, savepath=output_path("comparison", "combinedMultiplicityWithDomeCutPerGroupStats.png"),
    negLnLByGroup=negLnLByGroupWithDomeCut, chi2ByGroup=chi2ByGroupWithDomeCut,
    sourceLabel=r'$^{252}$Cf Source:' + '\nw/ dome cut', showPerGroupStats=True,
)

globalNormFactorWithoutDomeCut = global_normalization_factor(groupsWithoutDomeCut)
negLnLByGroupWithoutDomeCut = [
    neg_ln_l_calc(g, g["seitz"], globalNormFactorWithoutDomeCut, []) for g in groupsWithoutDomeCut
]
chi2_debug("[chi2 debug] -- without dome cut --")
chi2ByGroupWithoutDomeCut = [
    chi_squared_diagnostic(g, g["seitz"], globalNormFactorWithoutDomeCut, [], debug=True)
    for g in groupsWithoutDomeCut
]
plot_combined_multiplicity_single(
    groupsWithoutDomeCut, savepath=output_path("comparison", "combinedMultiplicityWithoutDomeCut.png"),
    negLnLByGroup=negLnLByGroupWithoutDomeCut, chi2ByGroup=chi2ByGroupWithoutDomeCut,
    sourceLabel=r'$^{252}$Cf Source:' + '\nw/o dome cut',
)
plot_combined_multiplicity_single(
    groupsWithoutDomeCut, savepath=output_path("comparison", "combinedMultiplicityWithoutDomeCutPerGroupStats.png"),
    negLnLByGroup=negLnLByGroupWithoutDomeCut, chi2ByGroup=chi2ByGroupWithoutDomeCut,
    sourceLabel=r'$^{252}$Cf Source:' + '\nw/o dome cut', showPerGroupStats=True,
)

globalNormFactorNoCutAtAll = global_normalization_factor(groupsNoCutAtAll)
negLnLByGroupNoCutAtAll = [
    neg_ln_l_calc(g, g["seitz"], globalNormFactorNoCutAtAll, []) for g in groupsNoCutAtAll
]
chi2_debug("[chi2 debug] -- no cut at all --")
chi2ByGroupNoCutAtAll = [
    chi_squared_diagnostic(g, g["seitz"], globalNormFactorNoCutAtAll, [], debug=True)
    for g in groupsNoCutAtAll
]
plot_combined_multiplicity_single(
    groupsNoCutAtAll, savepath=output_path("comparison", "combinedMultiplicityNoCutAtAll.png"),
    negLnLByGroup=negLnLByGroupNoCutAtAll, chi2ByGroup=chi2ByGroupNoCutAtAll,
    sourceLabel=r'$^{252}$Cf Source:' + '\nno cut at all',
)
plot_combined_multiplicity_single(
    groupsNoCutAtAll, savepath=output_path("comparison", "combinedMultiplicityNoCutAtAllPerGroupStats.png"),
    negLnLByGroup=negLnLByGroupNoCutAtAll, chi2ByGroup=chi2ByGroupNoCutAtAll,
    sourceLabel=r'$^{252}$Cf Source:' + '\nno cut at all', showPerGroupStats=True,
)

"""
SINGLE THRESHOLD RATES WITH AND WITHOUT THEORETICAL THRESHOLD
"""
binLabels3 = ["1", "2", "3+"]
binLabelsFull = ["1", "2", "3", "4", "5+"]

def plot_linhist(binLabels, binCounts, binCountError, backBins, backErrorLow, backErrorHigh,
                  backSubBins, backSubErrorLow, backSubErrorHigh, zeroKevRate, seitzRate, seitz,
                  savepath, bestFitRate=None, bestThreshold=None, bestNegLnL=None):
    plt.figure(figsize=(10, 10))
    x = np.arange(len(binLabels))
    edges = np.concatenate(([x[0] - 0.5], (x[:-1] + x[1:]) / 2, [x[-1] + 0.5]))

    plt.errorbar(x, binCounts, yerr=binCountError, fmt='o', color="red", ecolor="red", label="Source Rate")
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
        plot_z_distribution(
            singleBubbleZByGroup.get((g["p"], g["T"]), []), backgroundZByGroup.get((g["p"], g["T"]), []), g["seitz"],
            savepath=output_path(outputDir, f"zDistributions/zdist{g['p']}{g['T']}.png"),
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
            savepath=output_path(outputDir, f"negLnLScans/negLnLScan{g['p']}{g['T']}.png"),
        )
        if compareNormalizationModes:
            plot_neg_ln_l_scan_comparison(
                g["bestFit"]["thresholdRange"], g["bestFit"]["negLnLCurve"],
                g["bestFit"]["threshold"], g["bestFit"]["negLnL"], normalization_mode_label(useGlobalNormalization),
                g["bestFit"]["altNegLnLCurve"], g["bestFit"]["altThreshold"], g["bestFit"]["altNegLnL"],
                normalization_mode_label(not useGlobalNormalization),
                savepath=output_path(outputDir, f"negLnLScans/negLnLScanCompare{g['p']}{g['T']}.png"),
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
            savepath=output_path(outputDir, f"negLnLScansNoSingles/negLnLScanNoSingles{g['p']}{g['T']}.png"),
        )
        plot_neg_ln_l_scan_comparison(
            g["bestFit"]["thresholdRange"], g["bestFit"]["negLnLCurve"],
            g["bestFit"]["threshold"], g["bestFit"]["negLnL"], "All bins (1, 2, 3+)",
            g["bestFitNoSingles"]["negLnLCurve"], g["bestFitNoSingles"]["threshold"], g["bestFitNoSingles"]["negLnL"],
            "No singles (2, 3+ only)",
            savepath=output_path(outputDir, f"negLnLScansNoSingles/negLnLScanCompare{g['p']}{g['T']}.png"),
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
            binLabels3, rebin(g["binCounts"]), rebin_errors(g["binCountError"]),
            rebin(g["backBins"]), rebin_errors(g["backErrorLow"]), rebin_errors(g["backErrorHigh"]),
            g["backSub"], g["errLow"], g["errHigh"],
            zeroKevRateRebinned, g["seitzRate"], g["seitz"],
        )
    seitzRateFull = [normalizationFactor * c for c in g["seitzCountsFull"]]
    return (
        binLabelsFull, g["binCounts"], g["binCountError"], g["backBins"], g["backErrorLow"], g["backErrorHigh"],
        g["backSubFull"], g["backSubErrorLowFull"], g["backSubErrorHighFull"],
        zeroKevRateFull, seitzRateFull, g["seitz"],
    )


def plot_group_linhists(groups, zeroKevRateRebinned, zeroKevRateFull, normalizationFactor, outputDir):
    for g in groups:
        linhistArgs = build_linhist_args(g, zeroKevRateRebinned, zeroKevRateFull, normalizationFactor)
        bestFitRate = g["bestFit"]["rate"] if useRebinnedThresholdPlots else g["bestFit"]["rateFull"]

        # one version without a best-fit curve, one with the normal (all-bins) fit, one with the no-singles (2, 3+ only) fit
        plot_linhist(*linhistArgs, savepath=output_path(outputDir, f"linHists/linhist{g['p']}{g['T']}.png"))
        plot_linhist(
            *linhistArgs, savepath=output_path(outputDir, f"linHists/linhistFit{g['p']}{g['T']}.png"),
            bestFitRate=bestFitRate, bestThreshold=g["bestFit"]["threshold"], bestNegLnL=g["bestFit"]["negLnL"],
        )
        if computeNoSinglesFit:
            bestFitRateNoSingles = g["bestFitNoSingles"]["rate"] if useRebinnedThresholdPlots else g["bestFitNoSingles"]["rateFull"]
            plot_linhist(
                *linhistArgs, savepath=output_path(outputDir, f"linHists/linhistFitNoSingles{g['p']}{g['T']}.png"),
                bestFitRate=bestFitRateNoSingles, bestThreshold=g["bestFitNoSingles"]["threshold"],
                bestNegLnL=g["bestFitNoSingles"]["negLnL"],
            )


# the "average across every (p,T) group" version of the pipeline above: one combined rate group, one normal fit, one no-singles fit, three linhist variants
def run_avg_group_pipeline(groups, globalNormFactor, simExcludedRegions, positionDomeFlags,
                            zeroKevRateRebinned, zeroKevRateFull, outputDir):
    is_region_excluded = make_is_region_excluded(positionDomeFlags)
    # the avg group pools every (p, T)'s source events, so it pools their backgrounds too
    avgRateGroup = compute_rate_group(keep=lambda i, region: not is_region_excluded(i, region),
                                      backgroundPTs=pToUse)
    avgGroup = {
        **avgRateGroup,
        "backSub": rebin(avgRateGroup["backSubFull"]),
        "errLow": rebin_errors(avgRateGroup["backSubErrorLowFull"]),
        "errHigh": rebin_errors(avgRateGroup["backSubErrorHighFull"]),
        "backSubUnfloored": rebin(avgRateGroup["backSubUnflooredFull"]),
        "errLowUnfloored": rebin_errors(avgRateGroup["errLowUnflooredFull"]),
        "errHighUnfloored": rebin_errors(avgRateGroup["errHighUnflooredFull"]),
    }

    avgSeitz = np.mean([g["seitz"] for g in groups])
    avgSeitzCountsRaw, _ = get_multiplicity_counts(avgSeitz, excludedRegionsOverride=simExcludedRegions)

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
            binLabels3, rebin(avgGroup["binCounts"]), rebin_errors(avgGroup["binCountError"]),
            rebin(avgGroup["backBins"]), rebin_errors(avgGroup["backErrorLow"]), rebin_errors(avgGroup["backErrorHigh"]),
            avgGroup["backSub"], avgGroup["errLow"], avgGroup["errHigh"],
            zeroKevRateRebinned, avgSeitzRate, avgSeitz,
        )
        avgBestFitRate = avgBestFit["rate"]
    else:
        totalAvg = sum(avgGroup["backSubUnflooredFull"])
        avgSeitzRatiosFull = counts_to_ratios(avgSeitzCountsRaw)
        avgLinhistArgs = (
            binLabelsFull, avgGroup["binCounts"], avgGroup["binCountError"], avgGroup["backBins"],
            avgGroup["backErrorLow"], avgGroup["backErrorHigh"],
            avgGroup["backSubFull"], avgGroup["backSubErrorLowFull"], avgGroup["backSubErrorHighFull"],
            zeroKevRateFull, seitz_count(avgSeitzRatiosFull, totalAvg), avgSeitz,
        )
        avgBestFitRate = avgBestFit["rateFull"]

    # one version without a best-fit curve, one with the normal (all-bins) fit, one with the no-singles (2, 3+ only) fit
    plot_linhist(*avgLinhistArgs, savepath=output_path(outputDir, "avgseitz.png"))
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


def run_pipeline(groups, normalizationFactor, simExcludedRegions, positionDomeFlags, outputDir):
    plot_group_z_distributions(groups, outputDir)

    globalNormFactor = global_normalization_factor(groups)
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
                            zeroKevRateRebinned, zeroKevRateFull, outputDir)


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

    # 3 rebinned (1, 2, 3+) or 5 full (1, 2, 3, 4, 5+) multiplicity bins per (p, T) group
    binsPerGroup = 3 if useRebinnedThresholdPlots else 5
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
            savepath=output_path(outputDir, f"linHistsA/negLnLScanA{g['p']}{g['T']}.png"),
            bestLabel=f"A * Seitz fit: {g['bestFitA']['threshold']:0.2f} keV (" + r"$-2\Delta\ln L$" + f"={g['bestFitA']['negLnL']:0.2f})",
        )

        linhistArgs = build_linhist_args(g, zeroKevRateRebinned, zeroKevRateFull, normalizationFactor)
        bestFitRateA = g["bestFitA"]["rate"] if useRebinnedThresholdPlots else g["bestFitA"]["rateFull"]
        plot_linhist(
            *linhistArgs, savepath=output_path(outputDir, f"linHistsA/linhistFitA{g['p']}{g['T']}.png"),
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


run_pipeline(groupsWithoutDomeCut, normalizationFactorWithoutDomeCut, [], positionDomeFlagsPre, "withoutDomeCut")
run_pipeline(groupsWithDomeCut, normalizationFactorWithDomeCut, ["dome"], positionDomeFlagsPost, "withDomeCut")

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
        savepath=output_path("comparison", f"negLnLScans/negLnLScanDomeCutCompare{gWithCut['p']}{gWithCut['T']}.png"),
    )

# single-A "threshold = A * seitz" fit, run last: once for pre-cut (withoutDomeCut), once for post-cut (withDomeCut)
fit_seitz_multiplier(groupsWithoutDomeCut, [], "withoutDomeCut", normalizationFactorWithoutDomeCut)
fit_seitz_multiplier(groupsWithDomeCut, ["dome"], "withDomeCut", normalizationFactorWithDomeCut)
