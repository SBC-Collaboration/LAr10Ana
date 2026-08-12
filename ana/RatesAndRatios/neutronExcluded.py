## dome exclusion: same rate/threshold-fit pipeline as neutronPlotting.py, but the sim and real sides both drop dome-region hits, gated by excludedRegions below
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

# withoutDomeCut/withDomeCut pipeline (build_groups()/run_pipeline()): PRE is the always-applied loose baseline, POST adds the full dome cut on top -- independent of DOME_Z_THRESHOLD_CM above
DOME_Z_THRESHOLD_CM_PRE = 6
DOME_Z_THRESHOLD_CM_POST = -3


# everything lands under one "plots" root: plots/withoutDomeCut/, plots/withDomeCut/, plots/comparison/
PLOTS_ROOT = "plots"


def output_path(outputDir, filename):
    path = os.path.join(PLOTS_ROOT, outputDir, filename)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path


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
backgroundRunsWarm = ["20251113_9","20251113_10","20251113_11","20251114_0","20251114_1","20251114_6","20251114_36","20251114_37","20251115_0","20251115_1","20251115_2","20251115_3","20251115_4","20251115_5","20251116_1","20251116_2","20251117_0","20251117_1","20251126_7","20251126_8","20251127_0","20251127_1","20251127_2","20251127_3","20251127_4","20251127_5","20251128_0","20251128_1","20251128_2","20251128_3","20251128_4","20251129_0","20251129_1","20251129_2","20251129_3","20251129_4","20251129_5","20251130_0","20251130_1","20251130_2","20251130_3","20251130_4","20251130_5",]
#backgroundRunsWarm = []
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
neutronRunsWarmB = ["20260108_4", "20260108_5", "20260108_6", "20260108_7", "20260108_8", "20260109_0"]
#neutronRunsWarmB = []
## cold annular
neutronRunsColdB = ["20260122_3","20260122_4","20260122_5","20260122_6","20260123_0","20260123_1","20260123_2","20260123_3","20260123_4","20260123_8","20260123_9","20260123_10","20260124_0","20260124_1","20260124_3","20260124_4","20260124_5","20260125_0","20260125_1","20260125_2","20260125_3","20260125_4","20260125_5","20260125_6","20260125_7","20260125_8"]
## 119K
neutronRunsHotB = ["20260205_12","20260205_13","20260205_14","20260205_15","20260205_16","20260205_17","20260205_18","20260206_0","20260206_1","20260206_2","20260206_3","20260206_4","20260206_5","20260206_6","20260206_7","20260213_1","20260213_2","20260213_3","20260213_4","20260213_5","20260213_6","20260213_7","20260213_8","20260213_9","20260214_0","20260214_1","20260214_2","20260214_3","20260214_4","20260214_5","20260214_6","20260214_7","20260214_8","20260214_9","20260214_10","20260214_11","20260214_12","20260214_13","20260214_14","20260215_0","20260215_1","20260215_2","20260215_3","20260215_4","20260215_5","20260215_6","20260215_7","20260215_8","20260215_9","20260215_10","20260215_11","20260215_12","20260215_13","20260215_14","20260216_0","20260216_1","20260216_2","20260216_3","20260216_4","20260216_5","20260216_6","20260216_7","20260216_8","20260216_9","20260216_10","20260216_11","20260216_12","20260216_13","20260217_0","20260217_1","20260217_2","20260217_3","20260217_4","20260217_5","20260217_6"]

## ones that are used for this graph
useConfigB = True

neutronRuns = (neutronRunsWarmB + neutronRunsColdB + neutronRunsHotB) if useConfigB \
    else (neutronRunsWarm + neutronRunsCold + neutronRunsHot)

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
    totalExcluded = regionExcludedCount + zExcludedCount
    verb = "excluded" if applyDomeExclusion else "flagged (not excluded, kept in output)"
    print(f"[{label}] {totalExcluded} single-bubble events {verb} "
          f"({regionExcludedCount} by region label, {zExcludedCount} by reco Z)")
    return positions


# total background bin counts (multiplicity 1,2,3,4,5+) and live time
def load_background():
    subdirs = {os.path.basename(p.rstrip(os.sep)) for p in glob.glob(os.path.join(RECON_DIR, '*/'))}
    binCounts = [0] * 5
    liveTime = 0.0
    for run, ev, mult, _ in iter_matched_events(HANDSCAN_DIR, backgroundList):
        if run not in subdirs:
            continue
        expData = Streamer(os.path.join(RECON_DIR, run, 'exposure.sbc')).to_dict()
        for i in range(len(expData["ev"])):
            if int(expData["ev"][i]) == int(ev) and float(expData['PT2121_livetime'][i]) > 1:
                if mult != 0:
                    binCounts[min(mult, 5) - 1] += 1
                liveTime += float(expData['PT2121_livetime'][i])
                break
    return binCounts, liveTime

# (mult, region) per neutron-run event, its live time, and its (pset_lo, pset_hi, temp)
def load_neutron_events():
    bubbleCount, sourceTimes, psetsTemps, runEvs = [], [], [], []
    for run, ev, mult, region in iter_matched_events(HANDSCAN_DIR, neutronRuns):
        bubbleCount.append((mult, region))
        runEvs.append((run, ev))
        expData = Streamer(f'{RECON_DIR}{run}/exposure.sbc').to_dict()
        for i in range(len(expData["ev"])):
            if int(expData["ev"][i]) == int(ev):
                sourceTimes.append(float(expData['PT2121_livetime'][i]))
                break
        evData = Streamer(f'{RECON_DIR}{run}/event.sbc').to_dict()
        for i in range(len(evData["ev"])):
            if int(evData["ev"][i]) == int(ev):
                temp = 119 if (run in neutronRunsHot) else 116
                psetsTemps.append((evData["pset_lo"][i], evData["pset_hi"][i], temp))
                break
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
    backBins = [b * sourceTime / backgroundTime for b in backgroundBinCounts]
    backErrorLow = [np.sqrt(b) if b >= 1 else 0.0 for b in backBins]
    backErrorHigh = [np.sqrt(b) if b >= 1 else -np.log(1 - 0.68) / backgroundTime for b in backBins]

    binCountError = [np.sqrt(c) for c in binCounts]
    backSubBins = [c - b for c, b in zip(binCounts, backBins)]
    # background high -> subtracted rate pulled down; background low -> subtracted rate pulled up
    backSubErrorLow = [np.sqrt(countErr**2 + backErr**2) for countErr, backErr in zip(binCountError, backErrorHigh)]
    backSubErrorHigh = [np.sqrt(countErr**2 + backErr**2) for countErr, backErr in zip(binCountError, backErrorLow)]
    return backBins, backErrorLow, backErrorHigh, backSubBins, backSubErrorLow, backSubErrorHigh, binCountError


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
backgroundBinCounts, backgroundTime = load_background()
bubbleCount, sourceTimes, psetsTemps, neutronRunEvs = load_neutron_events()
singleBubblePositions = load_single_bubble_positions(neutronRuns, label="source")

# reco-Z-based dome flags for the withoutDomeCut/withDomeCut pipeline: "pre" = always-applied loose threshold, "post" = tighter threshold on top
positionDomeFlagsPre = compute_position_dome_flags(bubbleCount, neutronRunEvs, DOME_Z_THRESHOLD_CM_PRE)
positionDomeFlagsPost = compute_position_dome_flags(bubbleCount, neutronRunEvs, DOME_Z_THRESHOLD_CM_POST)

# (p, T) pairs with a fixed pressure setpoint, deduped
pToUse = sorted({(float(lo), float(t)) for lo, hi, t in psetsTemps if float(lo) == float(hi)})


# bins events matching `keep`, background-subtracts, and converts to a rate in counts/minute -- shared core of both build_groups()'s per-(p,T) loop and the avg-group section
def compute_rate_group(keep):
    binCounts, liveTimeSec = bin_multiplicities(bubbleCount, sourceTimes, keep)
    backBins, backErrorLow, backErrorHigh, backSubBins, backSubErrorLow, backSubErrorHigh, binCountError = \
        background_subtract(binCounts, liveTimeSec, backgroundBinCounts, backgroundTime)

    liveTimeMin = liveTimeSec / 60
    toRate = lambda values: [v / liveTimeMin for v in values]
    return {
        "liveTime": liveTimeMin,
        "binCounts": toRate(binCounts),
        "binCountError": toRate(binCountError),
        "backBins": toRate(backBins),
        "backErrorLow": toRate(backErrorLow),
        "backErrorHigh": toRate(backErrorHigh),
        "backSubFull": toRate(backSubBins),
        "backSubErrorLowFull": toRate(backSubErrorLow),
        "backSubErrorHighFull": toRate(backSubErrorHigh),
    }


# builds the full per-(p,T) "groups" list for one (simExcludedRegions, positionDomeFlags) combination -- lets the pre/post pipeline independently control the sim-side cut and real-side position cut
def build_groups(simExcludedRegions, positionDomeFlags):
    groups = []
    is_region_excluded = make_is_region_excluded(positionDomeFlags)

    for p, T in pToUse:
        rateGroup = compute_rate_group(
            keep=lambda i, region, p=p, T=T: not is_region_excluded(i, region)
                                              and psetsTemps[i][0] == p and psetsTemps[i][2] == T
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
            "seitzRate": rebin(seitzCounts),
            "seitzCounts": rebin(seitzCounts),
            "bestFit": {},
            **rateGroup,
        })

    # sort from lowest to highest seitz threshold
    groups.sort(key=lambda g: g["seitz"])

    # scale sim-predicted counts to match observed data rates, averaged across every group's own (data total)/(sim total) ratio
    observedRatesByGroup = [g["backSub"] for g in groups]
    simCountsByGroup = [g["seitzCounts"] for g in groups]
    normalizationFactor = np.mean([
        sum(observed) / sum(sim) for observed, sim in zip(observedRatesByGroup, simCountsByGroup)
    ])
    for g in groups:
        g["seitzRate"] = [normalizationFactor * ratio for ratio in g["seitzCounts"]]

    return groups, normalizationFactor


# "withoutDomeCut": sim cut off, real side at the loose threshold. "withDomeCut": sim cut on, real side at the tighter threshold. Each gets its own run_pipeline() subfolder
groupsWithoutDomeCut, normalizationFactorWithoutDomeCut = build_groups([], positionDomeFlagsPre)
groupsWithDomeCut, normalizationFactorWithDomeCut = build_groups(["dome"], positionDomeFlagsPost)

# dumps every group's background-subtracted data rate and Seitz/sim-predicted rate, pre-cut vs post-cut, to a plain text file
def write_rates_txt(groupsPre, groupsPost, savepath):
    binLabels = ["1", "2", "3+"]
    with open(savepath, "w") as f:
        f.write(f"Rates pre-cut (uncut baseline: real side Z > {DOME_Z_THRESHOLD_CM_PRE}cm, sim off) "
                f"vs post-cut (real side Z > {DOME_Z_THRESHOLD_CM_POST}cm, sim on), count/min\n")
        f.write("data = background-subtracted real rate, sim = Seitz/sim-predicted rate\n")
        f.write("=" * 78 + "\n\n")
        for gPre, gPost in zip(groupsPre, groupsPost):
            f.write(f"Seitz = {gPre['seitz']:0.2f} keV  (P = {gPre['p']:0.2f} bar, T = {gPre['T']:0.0f} K)\n")
            f.write(f"  Pre-cut  (Z > {DOME_Z_THRESHOLD_CM_PRE}cm, sim off):  livetime = {gPre['liveTime']:0.2f} min\n")
            for label, rate, errLow, errHigh, simRate in zip(
                    binLabels, gPre["backSub"], gPre["errLow"], gPre["errHigh"], gPre["seitzRate"]):
                f.write(f"    mult {label:<3}: data = {rate:0.4f} (+{errHigh:0.4f}/-{errLow:0.4f})   sim = {simRate:0.4f}\n")
            f.write(f"  Post-cut (Z > {DOME_Z_THRESHOLD_CM_POST}cm, sim on): livetime = {gPost['liveTime']:0.2f} min\n")
            for label, rate, errLow, errHigh, simRate in zip(
                    binLabels, gPost["backSub"], gPost["errLow"], gPost["errHigh"], gPost["seitzRate"]):
                f.write(f"    mult {label:<3}: data = {rate:0.4f} (+{errHigh:0.4f}/-{errLow:0.4f})   sim = {simRate:0.4f}\n")
            f.write("\n")
    print(f"wrote pre/post-cut rates to {savepath}")

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

plot_combined_multiplicity_comparison(
    groupsWithDomeCut, groupsWithoutDomeCut,
    savepath=output_path("comparison", "combinedMultiplicityComparison.png"),
)

"""
1-BUBBLE Z DISTRIBUTION, FOUR LOWEST SEITZ THRESHOLDS
"""
# how many of the lowest-seitz groups to plot
numLowestSeitzZDist = 4

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
            temp = 119 if (run in neutronRunsHot) else 116
            return float(evData["pset_lo"][i]), float(evData["pset_hi"][i]), temp
    return None


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
    totalExcluded = regionExcludedCount + zExcludedCount
    verb = "excluded" if applyDomeExclusion else "flagged (not excluded, kept in output)"
    print(f"[z-distribution source] {totalExcluded} single-bubble events {verb} "
          f"({regionExcludedCount} by region label, {zExcludedCount} by reco Z)")
    return byGroup


# z values below this are non-physical, so leave them out of the distribution plots
MIN_PHYSICAL_Z_MM = -300


# source/background counts per Z-bin, each divided by its own livetime [minutes] -> count/min
def plot_z_distribution(sourceZ, backgroundZ, seitz, savepath, sourceLiveTime, backgroundLiveTime, numBins=50 // 4):
    sourceZ = [z for z in sourceZ if z >= MIN_PHYSICAL_Z_MM]
    backgroundZ = [z for z in backgroundZ if z >= MIN_PHYSICAL_Z_MM]
    allZ = sourceZ + backgroundZ
    edges = np.linspace(min(allZ), max(allZ), numBins + 1) if allZ else numBins

    backgroundCounts, _ = np.histogram(backgroundZ, bins=edges)
    sourceCounts, _ = np.histogram(sourceZ, bins=edges)
    backgroundRate = backgroundCounts / backgroundLiveTime
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
backgroundSingleBubbleZ = [pos[4] for pos in load_single_bubble_positions(backgroundList, label="background", applyDomeExclusion=False)]
backgroundLiveTimeMin = backgroundTime / 60

# if true, use 1,2,3+ if false use 1,2,3,4,5+
useRebinnedThresholdPlots = True

"""
THEORETICAL THRESHOLDS
"""
# range to check for ratio matching
thresholdRange = np.arange(0.2, 30.1, 0.1)

# range of multipliers to scan for the single-A "threshold = A * seitz" fit below
seitzMultiplierRange = np.arange(0.1, 3.01, 0.01)

# chi-squared normalization mode: if True, scale every group's predicted counts by one shared factor across the whole dataset; if False, scale each to its own observed total (like oldCode.py)
useGlobalChi2Normalization = False

# if True, also compute the *other* normalization mode and generate the *Compare* plots overlaying both. If False (default), just use useGlobalChi2Normalization above
compareNormalizationModes = False

# if True, also fit/plot the no-singles (2, 3+ only) variant alongside the normal fit. Off by default since it's almost never needed
computeNoSinglesFit = False

def group_observed_total(dataGroup):
    rate = dataGroup["backSub"] if useRebinnedThresholdPlots else dataGroup["backSubFull"]
    return sum(v * dataGroup["liveTime"] for v in rate)

# livetime/total-weighted average of each group's own normalization, relative to the mean -- pulls a group's predicted total toward the dataset-wide average instead of forcing an exact match
def global_normalization_factor(dataGroups):
    weights = [group_observed_total(g) for g in dataGroups]
    meanWeight = np.mean(weights)
    avgRatio = sum(w ** 2 for w in weights) / sum(weights)
    return avgRatio / meanWeight

def normalization_mode_label(useGlobalNorm):
    return "Global normalization" if useGlobalNorm else "Per-threshold normalization"

# globalChi2NormFactor and simExcludedRegions are per-side, since they depend on which groups list / sim cut is active
def chi_squared_calc(dataGroup, estThreshold, globalChi2NormFactor, simExcludedRegions, useGlobalNorm=None):
    if useGlobalNorm is None:
        useGlobalNorm = useGlobalChi2Normalization

    if useRebinnedThresholdPlots:
        rate, errLowRate, errHighRate = dataGroup["backSub"], dataGroup["errLow"], dataGroup["errHigh"]
        predictedCounts = rebin(get_multiplicity_counts(estThreshold, excludedRegionsOverride=simExcludedRegions)[0])
    else:
        rate = dataGroup["backSubFull"]
        errLowRate, errHighRate = dataGroup["backSubErrorLowFull"], dataGroup["backSubErrorHighFull"]
        predictedCounts = get_multiplicity_counts(estThreshold, excludedRegionsOverride=simExcludedRegions)[0]

    # make sure to convert everything to a count properly before matching
    liveTime = dataGroup["liveTime"]
    observed = [v * liveTime for v in rate]
    errLow = [v * liveTime for v in errLowRate]
    errHigh = [v * liveTime for v in errHighRate]

    targetTotal = globalChi2NormFactor * sum(observed) if useGlobalNorm else sum(observed)
    predicted = seitz_count(counts_to_ratios(predictedCounts), targetTotal)

    chi2 = 0.0
    for obs, eLow, eHigh, pred in zip(observed, errLow, errHigh, predicted):
        err = eHigh if obs >= pred else eLow
        if err == 0:
            continue
        chi2 += ((obs - pred) / err) ** 2
    return chi2

# same idea as chi_squared_calc, but drops the multiplicity==1 bin and fits only 2 and 3+ -- always rebinned, results go in g["bestFitNoSingles"], never g["bestFit"]
def chi_squared_calc_no_singles(dataGroup, estThreshold, simExcludedRegions):
    rate, errLowRate, errHighRate = dataGroup["backSub"], dataGroup["errLow"], dataGroup["errHigh"]
    predictedCounts = rebin(get_multiplicity_counts(estThreshold, excludedRegionsOverride=simExcludedRegions)[0])

    liveTime = dataGroup["liveTime"]
    observed = [v * liveTime for v in rate][1:]
    errLow = [v * liveTime for v in errLowRate][1:]
    errHigh = [v * liveTime for v in errHighRate][1:]
    predictedCounts = predictedCounts[1:]

    targetTotal = sum(observed)
    predicted = seitz_count(counts_to_ratios(predictedCounts), targetTotal)

    chi2 = 0.0
    for obs, eLow, eHigh, pred in zip(observed, errLow, errHigh, predicted):
        err = eHigh if obs >= pred else eLow
        if err == 0:
            continue
        chi2 += ((obs - pred) / err) ** 2
    return chi2


def chi2_confidence_interval(gridKev, chi2Curve, bestIdx):
    """threshold values where chi2(x) - min(chi2) crosses 1, on each side of the best fit"""
    minChi2 = chi2Curve[bestIdx]
    target = minChi2 + 1.0

    lowThreshold = gridKev[0]
    for i in range(bestIdx, 0, -1):
        if chi2Curve[i - 1] >= target:
            x0, x1 = gridKev[i - 1], gridKev[i]
            y0, y1 = chi2Curve[i - 1], chi2Curve[i]
            frac = (target - y0) / (y1 - y0) if y1 != y0 else 0.0
            lowThreshold = x0 + frac * (x1 - x0)
            break

    highThreshold = gridKev[-1]
    for i in range(bestIdx, len(chi2Curve) - 1):
        if chi2Curve[i + 1] >= target:
            x0, x1 = gridKev[i], gridKev[i + 1]
            y0, y1 = chi2Curve[i], chi2Curve[i + 1]
            frac = (target - y0) / (y1 - y0) if y1 != y0 else 0.0
            highThreshold = x0 + frac * (x1 - x0)
            break

    return lowThreshold, highThreshold


# given a chi2 value per entry of thresholdRange, locates the best-fit threshold, its 1-sigma bounds, and the best-fit rate curve -- shared by every fit in the pipeline
# builds the standard {threshold, thresholdErrLow, thresholdErrHigh, chi2, rate, rateFull} fit dict for an already-known threshold/chi2 -- rate is scaled to match dataGroup's own observed total, same as the 0keV/avg-seitz reference lines. Shared by compute_best_fit() below and fit_seitz_multiplier()'s per-group bestFitA
def build_fit_result(dataGroup, threshold, chi2, simExcludedRegions, thresholdErrLow=0.0, thresholdErrHigh=0.0):
    bestFitRatios = counts_to_ratios(
        get_multiplicity_counts(threshold, excludedRegionsOverride=simExcludedRegions)[0]
    )
    bestFitRateFull = seitz_count(bestFitRatios, sum(dataGroup["backSubFull"]))
    return {
        "threshold": threshold,
        "thresholdErrLow": thresholdErrLow,
        "thresholdErrHigh": thresholdErrHigh,
        "chi2": chi2,
        "rate": rebin(bestFitRateFull),
        "rateFull": bestFitRateFull,
    }


def compute_best_fit(dataGroup, chi2Curve, simExcludedRegions):
    bestIdx = int(np.argmin(chi2Curve))
    bestThreshold = thresholdRange[bestIdx]
    lowThreshold, highThreshold = chi2_confidence_interval(thresholdRange, chi2Curve, bestIdx)

    fit = build_fit_result(
        dataGroup, bestThreshold, chi2Curve[bestIdx], simExcludedRegions,
        thresholdErrLow=bestThreshold - lowThreshold, thresholdErrHigh=highThreshold - bestThreshold,
    )
    fit["thresholdRange"] = thresholdRange
    fit["chi2Curve"] = chi2Curve
    return fit

# chi-squared vs threshold scan
def plot_chi2_scan(gridKev, chi2Curve, bestThreshold, bestChi2, savepath, xlabel="Threshold [keV]", bestLabel=None):
    if bestLabel is None:
        bestLabel = f"Best fit: {bestThreshold:0.2f} keV (" + r"$\chi^2$" + f"={bestChi2:0.2f})"
    plt.figure(figsize=(8, 6))
    plt.plot(gridKev, chi2Curve, 'o', markersize=3, color="steelblue")
    plt.axvline(bestThreshold, color='red', linestyle='--', label=bestLabel)
    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel(r"$\chi^2$", fontsize=16)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig(savepath)
    plt.close()

# generic two-series chi2 scan overlay -- used for normalization-mode, no-singles, and dome-cut comparisons alike. seriesA is blue, seriesB red
def plot_chi2_scan_comparison(gridKev, chi2CurveA, bestThresholdA, bestChi2A, labelA,
                               chi2CurveB, bestThresholdB, bestChi2B, labelB, savepath):
    plt.figure(figsize=(8, 6))
    plt.plot(gridKev, chi2CurveA, 'o', markersize=3, color="steelblue", label=labelA)
    plt.axvline(bestThresholdA, color='steelblue', linestyle='--',
                label=f"{labelA} best fit: {bestThresholdA:0.2f} keV (" + r"$\chi^2$" + f"={bestChi2A:0.2f})")
    plt.plot(gridKev, chi2CurveB, 'o', markersize=3, color="red", label=labelB)
    plt.axvline(bestThresholdB, color='red', linestyle='--',
                label=f"{labelB} best fit: {bestThresholdB:0.2f} keV (" + r"$\chi^2$" + f"={bestChi2B:0.2f})")
    plt.xlabel("Threshold [keV]", fontsize=16)
    plt.ylabel(r"$\chi^2$", fontsize=16)
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

# generic two-series fit-vs-seitz overlay, same idea as plot_chi2_scan_comparison() above
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
SINGLE THRESHOLD RATES WITH AND WITHOUT THEORETICAL THRESHOLD
"""
binLabels3 = ["1", "2", "3+"]
binLabelsFull = ["1", "2", "3", "4", "5+"]

def plot_linhist(binLabels, binCounts, binCountError, backBins, backErrorLow, backErrorHigh,
                  backSubBins, backSubErrorLow, backSubErrorHigh, zeroKevRate, seitzRate, seitz,
                  savepath, bestFitRate=None, bestThreshold=None, bestChi2=None):
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
                   label=f"Best Fit Threshold\n({bestThreshold:0.2f} keV, " + r"$\chi^2$" + f"={bestChi2:0.2f})")

    plt.xticks(x, binLabels, fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlabel("Bubble Multiplicity", fontsize=20)
    plt.ylabel("Rate [count/min]", fontsize=20)
    plt.legend(fontsize=20)
    plt.tight_layout()
    plt.savefig(savepath)
    plt.close()

"""
FULL PER-SIDE PIPELINE -- Z distributions, chi2 fits, linhists, avg group. One function per stage; run_pipeline() at the bottom calls them in order, once per side.
"""

def plot_group_z_distributions(groups, outputDir):
    for g in groups[:numLowestSeitzZDist]:
        plot_z_distribution(
            singleBubbleZByGroup.get((g["p"], g["T"]), []), backgroundSingleBubbleZ, g["seitz"],
            savepath=output_path(outputDir, f"zDistributions/zdist{g['p']}{g['T']}.png"),
            sourceLiveTime=g["liveTime"], backgroundLiveTime=backgroundLiveTimeMin,
        )


# fills every group's g["bestFit"] (normal, all-bins 1/2/3+ fit); if compareNormalizationModes is on, also fills g["bestFit"]["alt*"] for the side-by-side comparison plots only
def fit_groups(groups, globalChi2NormFactor, simExcludedRegions):
    for g in groups:
        chi2Curve = [chi_squared_calc(g, threshold, globalChi2NormFactor, simExcludedRegions)
                     for threshold in thresholdRange]
        g["bestFit"] = compute_best_fit(g, chi2Curve, simExcludedRegions)

        if compareNormalizationModes:
            altChi2Curve = [chi_squared_calc(g, threshold, globalChi2NormFactor, simExcludedRegions,
                                              useGlobalNorm=not useGlobalChi2Normalization)
                             for threshold in thresholdRange]
            altFit = compute_best_fit(g, altChi2Curve, simExcludedRegions)
            g["bestFit"].update({
                "altChi2Curve": altFit["chi2Curve"],
                "altThreshold": altFit["threshold"],
                "altThresholdErrLow": altFit["thresholdErrLow"],
                "altThresholdErrHigh": altFit["thresholdErrHigh"],
                "altChi2": altFit["chi2"],
            })


def plot_group_chi2_scans(groups, outputDir):
    for g in groups:
        plot_chi2_scan(
            g["bestFit"]["thresholdRange"], g["bestFit"]["chi2Curve"],
            g["bestFit"]["threshold"], g["bestFit"]["chi2"],
            savepath=output_path(outputDir, f"chi2Scans/chi2scan{g['p']}{g['T']}.png"),
        )
        if compareNormalizationModes:
            plot_chi2_scan_comparison(
                g["bestFit"]["thresholdRange"], g["bestFit"]["chi2Curve"],
                g["bestFit"]["threshold"], g["bestFit"]["chi2"], normalization_mode_label(useGlobalChi2Normalization),
                g["bestFit"]["altChi2Curve"], g["bestFit"]["altThreshold"], g["bestFit"]["altChi2"],
                normalization_mode_label(not useGlobalChi2Normalization),
                savepath=output_path(outputDir, f"chi2Scans/chi2scanCompare{g['p']}{g['T']}.png"),
            )


# separate fitting path, drops the multiplicity==1 bin -- stored in g["bestFitNoSingles"], never touches g["bestFit"] from fit_groups() above
def fit_groups_no_singles(groups, simExcludedRegions):
    for g in groups:
        chi2CurveNoSingles = [chi_squared_calc_no_singles(g, threshold, simExcludedRegions)
                               for threshold in thresholdRange]
        g["bestFitNoSingles"] = compute_best_fit(g, chi2CurveNoSingles, simExcludedRegions)


def plot_group_no_singles_comparisons(groups, outputDir):
    for g in groups:
        plot_chi2_scan(
            g["bestFitNoSingles"]["thresholdRange"], g["bestFitNoSingles"]["chi2Curve"],
            g["bestFitNoSingles"]["threshold"], g["bestFitNoSingles"]["chi2"],
            savepath=output_path(outputDir, f"chi2ScansNoSingles/chi2scanNoSingles{g['p']}{g['T']}.png"),
        )
        plot_chi2_scan_comparison(
            g["bestFit"]["thresholdRange"], g["bestFit"]["chi2Curve"],
            g["bestFit"]["threshold"], g["bestFit"]["chi2"], "All bins (1, 2, 3+)",
            g["bestFitNoSingles"]["chi2Curve"], g["bestFitNoSingles"]["threshold"], g["bestFitNoSingles"]["chi2"],
            "No singles (2, 3+ only)",
            savepath=output_path(outputDir, f"chi2ScansNoSingles/chi2scanCompare{g['p']}{g['T']}.png"),
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


# 0keV reference-line shape (no energy cut at all), rebinned and full versions -- each linhist scales this to its own group's total, since on its own it's tied to no particular normalization
def compute_zero_kev_reference(simExcludedRegions):
    zeroKevCountsRaw, _ = get_multiplicity_counts(0.0, excludedRegionsOverride=simExcludedRegions)
    zeroKevRatiosRebinned = counts_to_ratios(rebin(zeroKevCountsRaw))
    zeroKevRatiosFull = counts_to_ratios(zeroKevCountsRaw)
    return zeroKevRatiosRebinned, zeroKevRatiosFull


# builds the positional args plot_linhist() expects for one group, respecting useRebinnedThresholdPlots -- shared by plot_group_linhists() and fit_seitz_multiplier()'s linHistsA plots
def build_linhist_args(g, zeroKevRatiosRebinned, zeroKevRatiosFull):
    if useRebinnedThresholdPlots:
        # 0keV reference scaled to this group's own observed total, same as the full/non-rebinned branch below
        zeroKevRate3 = seitz_count(zeroKevRatiosRebinned, sum(g["backSub"]))
        return (
            binLabels3, rebin(g["binCounts"]), rebin_errors(g["binCountError"]),
            rebin(g["backBins"]), rebin_errors(g["backErrorLow"]), rebin_errors(g["backErrorHigh"]),
            g["backSub"], g["errLow"], g["errHigh"],
            zeroKevRate3, g["seitzRate"], g["seitz"],
        )
    total = sum(g["backSubFull"])
    return (
        binLabelsFull, g["binCounts"], g["binCountError"], g["backBins"], g["backErrorLow"], g["backErrorHigh"],
        g["backSubFull"], g["backSubErrorLowFull"], g["backSubErrorHighFull"],
        seitz_count(zeroKevRatiosFull, total), seitz_count(counts_to_ratios(g["seitzCountsFull"]), total), g["seitz"],
    )


def plot_group_linhists(groups, zeroKevRatiosRebinned, zeroKevRatiosFull, outputDir):
    for g in groups:
        linhistArgs = build_linhist_args(g, zeroKevRatiosRebinned, zeroKevRatiosFull)
        bestFitRate = g["bestFit"]["rate"] if useRebinnedThresholdPlots else g["bestFit"]["rateFull"]

        # one version without a best-fit curve, one with the normal (all-bins) fit, one with the no-singles (2, 3+ only) fit
        plot_linhist(*linhistArgs, savepath=output_path(outputDir, f"linHists/linhist{g['p']}{g['T']}.png"))
        plot_linhist(
            *linhistArgs, savepath=output_path(outputDir, f"linHists/linhistFit{g['p']}{g['T']}.png"),
            bestFitRate=bestFitRate, bestThreshold=g["bestFit"]["threshold"], bestChi2=g["bestFit"]["chi2"],
        )
        if computeNoSinglesFit:
            bestFitRateNoSingles = g["bestFitNoSingles"]["rate"] if useRebinnedThresholdPlots else g["bestFitNoSingles"]["rateFull"]
            plot_linhist(
                *linhistArgs, savepath=output_path(outputDir, f"linHists/linhistFitNoSingles{g['p']}{g['T']}.png"),
                bestFitRate=bestFitRateNoSingles, bestThreshold=g["bestFitNoSingles"]["threshold"],
                bestChi2=g["bestFitNoSingles"]["chi2"],
            )


# the "average across every (p,T) group" version of the pipeline above: one combined rate group, one normal fit, one no-singles fit, three linhist variants
def run_avg_group_pipeline(groups, globalChi2NormFactor, simExcludedRegions, positionDomeFlags,
                            zeroKevRatiosRebinned, zeroKevRatiosFull, outputDir):
    is_region_excluded = make_is_region_excluded(positionDomeFlags)
    avgRateGroup = compute_rate_group(keep=lambda i, region: not is_region_excluded(i, region))
    avgGroup = {
        **avgRateGroup,
        "backSub": rebin(avgRateGroup["backSubFull"]),
        "errLow": rebin_errors(avgRateGroup["backSubErrorLowFull"]),
        "errHigh": rebin_errors(avgRateGroup["backSubErrorHighFull"]),
    }

    avgSeitz = np.mean([g["seitz"] for g in groups])
    avgSeitzCountsRaw, _ = get_multiplicity_counts(avgSeitz, excludedRegionsOverride=simExcludedRegions)

    avgChi2Curve = [chi_squared_calc(avgGroup, threshold, globalChi2NormFactor, simExcludedRegions)
                     for threshold in thresholdRange]
    avgBestFit = compute_best_fit(avgGroup, avgChi2Curve, simExcludedRegions)

    if compareNormalizationModes:
        avgAltChi2Curve = [chi_squared_calc(avgGroup, threshold, globalChi2NormFactor, simExcludedRegions,
                                             useGlobalNorm=not useGlobalChi2Normalization)
                            for threshold in thresholdRange]
        avgAltFit = compute_best_fit(avgGroup, avgAltChi2Curve, simExcludedRegions)
        avgBestFit.update({
            "altChi2Curve": avgAltFit["chi2Curve"],
            "altThreshold": avgAltFit["threshold"],
            "altChi2": avgAltFit["chi2"],
        })

    if computeNoSinglesFit:
        # separate no-singles fit for the avg group, same as fit_groups_no_singles() above
        avgChi2CurveNoSingles = [chi_squared_calc_no_singles(avgGroup, threshold, simExcludedRegions)
                                  for threshold in thresholdRange]
        avgBestFitNoSingles = compute_best_fit(avgGroup, avgChi2CurveNoSingles, simExcludedRegions)

        plot_chi2_scan(
            avgBestFitNoSingles["thresholdRange"], avgBestFitNoSingles["chi2Curve"],
            avgBestFitNoSingles["threshold"], avgBestFitNoSingles["chi2"],
            savepath=output_path(outputDir, "chi2scanAvgNoSingles.png"),
        )
        plot_chi2_scan_comparison(
            avgBestFit["thresholdRange"], avgBestFit["chi2Curve"],
            avgBestFit["threshold"], avgBestFit["chi2"], "All bins (1, 2, 3+)",
            avgBestFitNoSingles["chi2Curve"], avgBestFitNoSingles["threshold"], avgBestFitNoSingles["chi2"],
            "No singles (2, 3+ only)",
            savepath=output_path(outputDir, "chi2scanAvgCompareSingles.png"),
        )

    plot_chi2_scan(
        avgBestFit["thresholdRange"], avgBestFit["chi2Curve"],
        avgBestFit["threshold"], avgBestFit["chi2"],
        savepath=output_path(outputDir, "chi2scanAvg.png"),
    )
    if compareNormalizationModes:
        plot_chi2_scan_comparison(
            avgBestFit["thresholdRange"], avgBestFit["chi2Curve"],
            avgBestFit["threshold"], avgBestFit["chi2"], normalization_mode_label(useGlobalChi2Normalization),
            avgBestFit["altChi2Curve"], avgBestFit["altThreshold"], avgBestFit["altChi2"],
            normalization_mode_label(not useGlobalChi2Normalization),
            savepath=output_path(outputDir, "chi2scanCompareAvg.png"),
        )

    if useRebinnedThresholdPlots:
        avgSeitzCountsRebinned = rebin(avgSeitzCountsRaw)
        # same per-total scaling as plot_group_linhists() above -- normalizationFactor isn't reliable here since avgSeitz usually falls outside the thresholds it was calibrated on
        avgSeitzRate = seitz_count(counts_to_ratios(avgSeitzCountsRebinned), sum(avgGroup["backSub"]))
        avgZeroKevRate3 = seitz_count(zeroKevRatiosRebinned, sum(avgGroup["backSub"]))
        avgLinhistArgs = (
            binLabels3, rebin(avgGroup["binCounts"]), rebin_errors(avgGroup["binCountError"]),
            rebin(avgGroup["backBins"]), rebin_errors(avgGroup["backErrorLow"]), rebin_errors(avgGroup["backErrorHigh"]),
            avgGroup["backSub"], avgGroup["errLow"], avgGroup["errHigh"],
            avgZeroKevRate3, avgSeitzRate, avgSeitz,
        )
        avgBestFitRate = avgBestFit["rate"]
    else:
        totalAvg = sum(avgGroup["backSubFull"])
        avgSeitzRatiosFull = counts_to_ratios(avgSeitzCountsRaw)
        avgLinhistArgs = (
            binLabelsFull, avgGroup["binCounts"], avgGroup["binCountError"], avgGroup["backBins"],
            avgGroup["backErrorLow"], avgGroup["backErrorHigh"],
            avgGroup["backSubFull"], avgGroup["backSubErrorLowFull"], avgGroup["backSubErrorHighFull"],
            seitz_count(zeroKevRatiosFull, totalAvg), seitz_count(avgSeitzRatiosFull, totalAvg), avgSeitz,
        )
        avgBestFitRate = avgBestFit["rateFull"]

    # one version without a best-fit curve, one with the normal (all-bins) fit, one with the no-singles (2, 3+ only) fit
    plot_linhist(*avgLinhistArgs, savepath=output_path(outputDir, "avgseitz.png"))
    plot_linhist(
        *avgLinhistArgs, savepath=output_path(outputDir, "avgseitzFit.png"),
        bestFitRate=avgBestFitRate, bestThreshold=avgBestFit["threshold"], bestChi2=avgBestFit["chi2"],
    )
    if computeNoSinglesFit:
        avgBestFitRateNoSingles = avgBestFitNoSingles["rate"] if useRebinnedThresholdPlots else avgBestFitNoSingles["rateFull"]
        plot_linhist(
            *avgLinhistArgs, savepath=output_path(outputDir, "avgseitzFitNoSingles.png"),
            bestFitRate=avgBestFitRateNoSingles, bestThreshold=avgBestFitNoSingles["threshold"],
            bestChi2=avgBestFitNoSingles["chi2"],
        )


def run_pipeline(groups, normalizationFactor, simExcludedRegions, positionDomeFlags, outputDir):
    plot_group_z_distributions(groups, outputDir)

    globalChi2NormFactor = global_normalization_factor(groups)
    fit_groups(groups, globalChi2NormFactor, simExcludedRegions)
    plot_group_chi2_scans(groups, outputDir)

    plot_fit_to_seitz_ratio(groups, savepath=output_path(outputDir, "fitToSeitzRatio.png"))
    if compareNormalizationModes:
        plot_fit_to_seitz_ratio_comparison(
            [g["seitz"] for g in groups],
            *fit_threshold_series(groups, "bestFit"), normalization_mode_label(useGlobalChi2Normalization),
            *fit_threshold_series(groups, "bestFit", "altThreshold", "altThresholdErrLow", "altThresholdErrHigh"),
            normalization_mode_label(not useGlobalChi2Normalization),
            savepath=output_path(outputDir, "fitToSeitzRatioCompare.png"),
        )

    if computeNoSinglesFit:
        fit_groups_no_singles(groups, simExcludedRegions)
        plot_group_no_singles_comparisons(groups, outputDir)

    zeroKevRatiosRebinned, zeroKevRatiosFull = compute_zero_kev_reference(simExcludedRegions)
    plot_group_linhists(groups, zeroKevRatiosRebinned, zeroKevRatiosFull, outputDir)

    run_avg_group_pipeline(groups, globalChi2NormFactor, simExcludedRegions, positionDomeFlags,
                            zeroKevRatiosRebinned, zeroKevRatiosFull, outputDir)


# total chi2 across every group when every group is fit against threshold = A * g["seitz"], for one shared A instead of a separate best-fit threshold per group
def chi_squared_calc_seitz_multiplier(groups, A, globalChi2NormFactor, simExcludedRegions):
    return sum(chi_squared_calc(g, A * g["seitz"], globalChi2NormFactor, simExcludedRegions) for g in groups)


# single-parameter fit: scans A over seitzMultiplierRange, picks the A minimizing the combined chi2 across every (p,T) group at once, then reports it per group as threshold = A * g["seitz"] (stored in g["bestFitA"], alongside g["bestFit"]/g["bestFitNoSingles"]). Run once per side, after run_pipeline(). All plots land in outputDir/linHistsA/
def fit_seitz_multiplier(groups, simExcludedRegions, outputDir):
    globalChi2NormFactor = global_normalization_factor(groups)
    chi2Curve = [chi_squared_calc_seitz_multiplier(groups, A, globalChi2NormFactor, simExcludedRegions)
                 for A in seitzMultiplierRange]
    bestIdx = int(np.argmin(chi2Curve))
    bestA = seitzMultiplierRange[bestIdx]
    lowA, highA = chi2_confidence_interval(seitzMultiplierRange, chi2Curve, bestIdx)
    bestChi2 = chi2Curve[bestIdx]
    AErrLow, AErrHigh = bestA - lowA, highA - bestA

    plot_chi2_scan(
        seitzMultiplierRange, chi2Curve, bestA, bestChi2,
        savepath=output_path(outputDir, "chi2scanSeitzMultiplier.png"),
        xlabel="Seitz threshold multiplier A",
        bestLabel=f"Best fit: A = {bestA:0.3f} (" + r"$\chi^2$" + f"={bestChi2:0.2f})",
    )
    print(f"[{outputDir}] best-fit Seitz threshold multiplier: A = {bestA:0.3f} "
          f"(+{AErrHigh:0.3f}/-{AErrLow:0.3f}), chi2 = {bestChi2:0.2f}")

    # per-group threshold = A * seitz, with its own chi2 and A's uncertainty propagated multiplicatively
    for g in groups:
        threshold = bestA * g["seitz"]
        groupChi2 = chi_squared_calc(g, threshold, globalChi2NormFactor, simExcludedRegions)
        g["bestFitA"] = build_fit_result(
            g, threshold, groupChi2, simExcludedRegions,
            thresholdErrLow=AErrLow * g["seitz"], thresholdErrHigh=AErrHigh * g["seitz"],
        )

    # same trio of plots the normal fit gets (per-group chi2 scan, fit-vs-seitz scatter, linhist with the fit overlaid), all under linHistsA/
    plot_fit_to_seitz_ratio(groups, savepath=output_path(outputDir, "linHistsA/fitToSeitzRatioA.png"), bestFitKey="bestFitA")

    zeroKevRatiosRebinned, zeroKevRatiosFull = compute_zero_kev_reference(simExcludedRegions)
    for g in groups:
        plot_chi2_scan(
            g["bestFit"]["thresholdRange"], g["bestFit"]["chi2Curve"],
            g["bestFitA"]["threshold"], g["bestFitA"]["chi2"],
            savepath=output_path(outputDir, f"linHistsA/chi2scanA{g['p']}{g['T']}.png"),
            bestLabel=f"A * Seitz fit: {g['bestFitA']['threshold']:0.2f} keV (" + r"$\chi^2$" + f"={g['bestFitA']['chi2']:0.2f})",
        )

        linhistArgs = build_linhist_args(g, zeroKevRatiosRebinned, zeroKevRatiosFull)
        bestFitRateA = g["bestFitA"]["rate"] if useRebinnedThresholdPlots else g["bestFitA"]["rateFull"]
        plot_linhist(
            *linhistArgs, savepath=output_path(outputDir, f"linHistsA/linhistFitA{g['p']}{g['T']}.png"),
            bestFitRate=bestFitRateA, bestThreshold=g["bestFitA"]["threshold"], bestChi2=g["bestFitA"]["chi2"],
        )

    return {
        "aRange": seitzMultiplierRange,
        "chi2Curve": chi2Curve,
        "A": bestA,
        "AErrLow": AErrLow,
        "AErrHigh": AErrHigh,
        "chi2": bestChi2,
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

# same comparison, broken out per (p,T) group as chi2-vs-threshold scans
for gWithCut, gWithoutCut in zip(groupsWithDomeCut, groupsWithoutDomeCut):
    plot_chi2_scan_comparison(
        gWithCut["bestFit"]["thresholdRange"], gWithCut["bestFit"]["chi2Curve"],
        gWithCut["bestFit"]["threshold"], gWithCut["bestFit"]["chi2"], "Dome cut",
        gWithoutCut["bestFit"]["chi2Curve"], gWithoutCut["bestFit"]["threshold"], gWithoutCut["bestFit"]["chi2"],
        "No dome cut",
        savepath=output_path("comparison", f"chi2Scans/chi2scanDomeCutCompare{gWithCut['p']}{gWithCut['T']}.png"),
    )

# single-A "threshold = A * seitz" fit, run last: once for pre-cut (withoutDomeCut), once for post-cut (withDomeCut)
fit_seitz_multiplier(groupsWithoutDomeCut, [], "withoutDomeCut")
fit_seitz_multiplier(groupsWithDomeCut, ["dome"], "withDomeCut")
