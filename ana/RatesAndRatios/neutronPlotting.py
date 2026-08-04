## imports, self explanitory
import glob
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from sbcbinaryformat import Streamer
import SeitzModel as sm
## this only works from my home folder right now since i need the simulation output processed with a function in this directory
CF_SIM_DIR = "/nashome/o/ochiarin/Documents/neutronSim"
sys.path.insert(0, CF_SIM_DIR)
from cfconfBThresholds import get_multiplicity_counts


## config variables
PLOTS_DIR = "plots"
HANDSCAN_DIR = "/exp/e961/data/SBC-25-handscan/"
RECON_DIR = "/exp/e961/data/SBC-25-recon/v0.4.2/"


def plot_path(filename):
    path = os.path.join(PLOTS_DIR, filename)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path

# background rate calculation for subtraction
## warm annular
backgroundRunsWarm = ["20251113_9","20251113_10","20251113_11","20251114_0","20251114_1","20251114_6","20251114_36","20251114_37","20251115_0","20251115_1","20251115_2","20251115_3","20251115_4","20251115_5","20251116_1","20251116_2","20251117_0","20251117_1","20251126_7","20251126_8","20251127_0","20251127_1","20251127_2","20251127_3","20251127_4","20251127_5","20251128_0","20251128_1","20251128_2","20251128_3","20251128_4","20251129_0","20251129_1","20251129_2","20251129_3","20251129_4","20251129_5","20251130_0","20251130_1","20251130_2","20251130_3","20251130_4","20251130_5",]
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
## cold annular
neutronRunsColdB = ["20260122_3","20260122_4","20260122_5","20260122_6","20260123_0","20260123_1","20260123_2","20260123_3","20260123_4","20260123_8","20260123_9","20260123_10","20260124_0","20260124_1","20260124_3","20260124_4","20260124_5","20260125_0","20260125_1","20260125_2","20260125_3","20260125_4","20260125_5","20260125_6","20260125_7","20260125_8"]
## 119K
neutronRunsHotB = ["20260205_12","20260205_13","20260205_14","20260205_15","20260205_16","20260205_17","20260205_18","20260206_0","20260206_1","20260206_2","20260206_3","20260206_4","20260206_5","20260206_6","20260206_7","20260213_1","20260213_2","20260213_3","20260213_4","20260213_5","20260213_6","20260213_7","20260213_8","20260213_9","20260214_0","20260214_1","20260214_2","20260214_3","20260214_4","20260214_5","20260214_6","20260214_7","20260214_8","20260214_9","20260214_10","20260214_11","20260214_12","20260214_13","20260214_14","20260215_0","20260215_1","20260215_2","20260215_3","20260215_4","20260215_5","20260215_6","20260215_7","20260215_8","20260215_9","20260215_10","20260215_11","20260215_12","20260215_13","20260215_14","20260216_0","20260216_1","20260216_2","20260216_3","20260216_4","20260216_5","20260216_6","20260216_7","20260216_8","20260216_9","20260216_10","20260216_11","20260216_12","20260216_13","20260217_0","20260217_1","20260217_2","20260217_3","20260217_4","20260217_5","20260217_6"]

## ones that are used for this graph
useConfigB = True

neutronRuns = (neutronRunsWarmB + neutronRunsColdB + neutronRunsHotB) if useConfigB \
    else (neutronRunsWarm + neutronRunsCold + neutronRunsHot)

## idk if we need this but could help
excludedRegions = []

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
    bubbleCount, sourceTimes, psetsTemps = [], [], []
    for run, ev, mult, region in iter_matched_events(HANDSCAN_DIR, neutronRuns):
        bubbleCount.append((mult, region))
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
    return bubbleCount, sourceTimes, psetsTemps

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
    backSubErrorLow = [np.sqrt(ce**2 + be**2) for ce, be in zip(binCountError, backErrorHigh)]
    backSubErrorHigh = [np.sqrt(ce**2 + be**2) for ce, be in zip(binCountError, backErrorLow)]
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

# load in data
backgroundBinCounts, backgroundTime = load_background()
bubbleCount, sourceTimes, psetsTemps = load_neutron_events()

# (p, T) pairs with a fixed pressure setpoint, deduped
pToUse = sorted({(float(lo), float(t)) for lo, hi, t in psetsTemps if float(lo) == float(hi)})

# dictionaries for quick lookup when making plots
groups = []

for p, T in pToUse:
    binCounts, sourceTime = bin_multiplicities(
        bubbleCount, sourceTimes,
        keep=lambda i, region: region not in excludedRegions and psetsTemps[i][0] == p and psetsTemps[i][2] == T,
    )
    backBins, backErrorLow, backErrorHigh, backSubBins, backSubErrorLow, backSubErrorHigh, binCountError = background_subtract(
        binCounts, sourceTime, backgroundBinCounts, backgroundTime
    )

    # convert counts from source to a rate in counts/minute
    sourceTime /= 60
    binCounts = [v / sourceTime for v in binCounts]
    binCountError = [v / sourceTime for v in binCountError]
    backBins = [v / sourceTime for v in backBins]
    backErrorLow = [v / sourceTime for v in backErrorLow]
    backErrorHigh = [v / sourceTime for v in backErrorHigh]
    backSubBins = [v / sourceTime for v in backSubBins]
    backSubErrorLow = [v / sourceTime for v in backSubErrorLow]
    backSubErrorHigh = [v / sourceTime for v in backSubErrorHigh]

    # seitz threshold for this (P,T) pair, fed straight into the Cf sim to get counts
    seitz = sm.SeitzModel(p * 14.5038, -273.15 + T, 'argon').Q
    seitzCounts, _ = get_multiplicity_counts(seitz)
    groups.append({
        "p": p,
        "T": T,
        "seitz": seitz,
        "liveTime": sourceTime,
        "binCounts": binCounts,
        "binCountError": binCountError,
        "backBins": backBins,
        "backErrorLow": backErrorLow,
        "backErrorHigh": backErrorHigh,
        "backSubFull": backSubBins,
        "backSubErrorLowFull": backSubErrorLow,
        "backSubErrorHighFull": backSubErrorHigh,
        "seitzCountsFull": seitzCounts,
        "backSub": rebin(backSubBins),
        "errLow": rebin_errors(backSubErrorLow),
        "errHigh": rebin_errors(backSubErrorHigh),
        "seitzRate": rebin(seitzCounts),
        "seitzCounts": rebin(seitzCounts),
        "bestFit": {},
    })

# sort from lowest to highest seitz threshold
groups.sort(key=lambda g: g["seitz"])


# normalize everything
# data rate per threshold i per multiplicity j
r_ij = [g["backSub"] for g in groups]
# seitz counts per threshold i per multiplicity j
s_ij = [g["seitzCounts"] for g in groups]
normalizationFactor = np.mean([sum(r) / sum(s) for r, s in zip(r_ij, s_ij)])
for g in groups:
    g["seitzRate"] = [normalizationFactor * ratio for ratio in g["seitzCounts"]]

## plot making
"""
COMBINED PAPER PLOT
"""
def plot_combined_multiplicity(groups, savepath, groupsPerRow=4):
    binLabels = ["1", "2", "3+"]
    numBins = 3
    barWidth = 1.0
    gap = 0.6
    colWidthInches = 3.5

    nRows = int(np.ceil(len(groups) / groupsPerRow))
    fig, axes = plt.subplots(nRows, 1, figsize=(colWidthInches, 3.2 * nRows), squeeze=False)
    axes = axes[:, 0]

    globalMax = max(max(max(g["seitzRate"]), max(b + e for b, e in zip(g["backSub"], g["errHigh"])))
                     for g in groups)

    # plotting and formattingm, with help from claude
    for rowIdx, ax in enumerate(axes):
        trans = ax.get_xaxis_transform()
        rowGroups = groups[rowIdx * groupsPerRow:(rowIdx + 1) * groupsPerRow]

        pos = 0
        for gi, g in enumerate(rowGroups):
            xs = np.arange(pos, pos + numBins)

            ax.bar(xs, g["seitzRate"], width=barWidth, color="lightblue", edgecolor="steelblue", zorder=1)
            ax.errorbar(xs, g["backSub"], yerr=[g["errLow"], g["errHigh"]], fmt='o', color="red",
                        ecolor="red", zorder=2, markersize=3, elinewidth=1, capsize=2)


            for x, label in zip(xs, binLabels):
                ax.text(x, -0.05, label, transform=trans, ha='center', va='top', fontsize=12)
            # seitz threshold label
            center = (xs[0] + xs[-1]) / 2
            ax.text(center, 0.97, f'{g["seitz"]:0.2f}', transform=trans, ha='center', va='top', fontsize=10)

            # seperator 
            if gi < len(rowGroups) - 1:
                ax.axvline(pos + numBins + gap / 2 - 0.5, linestyle='--', linewidth=0.7,
                           color='gray', zorder=0)

            pos += numBins + gap

        ax.set_ylim(0, globalMax * 1.2)
        ax.set_xlim(-1, pos - gap)
        ax.set_xticks([])
        ax.tick_params(axis='y', labelsize=12)

    axes[0].set_title(r'$Q_{seitz}$ [keV]', loc='left', fontsize=16, pad=2)
    axes[0].set_title(r'$^{252}$Cf Source', loc='right', fontsize=16, pad=2)
    axes[nRows // 2].set_ylabel("Rate [count/min]", fontsize=16)
    axes[-1].set_xlabel("Bubble Multiplicity", fontsize=12, labelpad=28)
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.15)
    fig.savefig(savepath)
    plt.close(fig)

plot_combined_multiplicity(groups, savepath=plot_path("combinedMultiplicity.png"))

# if true, use 1,2,3+ if false use 1,2,3,4,5+
useRebinnedThresholdPlots = True

"""
THEORETICAL THRESHOLDS
"""
# range to check for ratio matching
thresholdRange = np.arange(0.05, 60.05, 0.05)

# chi-squared normalization mode: if True, every group's predicted counts are scaled
# by one shared normalization factor computed across the whole dataset; if False,
# each group's predicted counts are scaled to match its own observed total, like the
# old code (see oldCode.py's seitz_counts_per_threshold vs seitz_counts_entire_dataset
# for the same per-group vs whole-dataset dichotomy)
useGlobalChi2Normalization = False

def group_observed_total(dataGroup):
    rate = dataGroup["backSub"] if useRebinnedThresholdPlots else dataGroup["backSubFull"]
    return sum(v * dataGroup["liveTime"] for v in rate)

# livetime/total-weighted average of each group's own normalization, relative to the
# mean group total - applying this to a group's own total pulls its predicted total
# toward the dataset-wide average instead of forcing an exact match
def global_normalization_factor(dataGroups):
    weights = [group_observed_total(g) for g in dataGroups]
    meanWeight = np.mean(weights)
    avgRatio = sum(w ** 2 for w in weights) / sum(weights)
    return avgRatio / meanWeight

globalChi2NormFactor = global_normalization_factor(groups)

def normalization_mode_label(useGlobalNorm):
    return "Global normalization" if useGlobalNorm else "Per-threshold normalization"

# get the chi squared
def chi_squared_calc(dataGroup, estThreshold, useGlobalNorm=None):
    if useGlobalNorm is None:
        useGlobalNorm = useGlobalChi2Normalization

    if useRebinnedThresholdPlots:
        rate, errLowRate, errHighRate = dataGroup["backSub"], dataGroup["errLow"], dataGroup["errHigh"]
        predictedCounts = rebin(get_multiplicity_counts(estThreshold)[0])
    else:
        rate = dataGroup["backSubFull"]
        errLowRate, errHighRate = dataGroup["backSubErrorLowFull"], dataGroup["backSubErrorHighFull"]
        predictedCounts = get_multiplicity_counts(estThreshold)[0]

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

for g in groups:
    # all the chi 2 for matched thresholds
    chi2Curve = [chi_squared_calc(g, threshold) for threshold in thresholdRange]
    bestIdx = int(np.argmin(chi2Curve))
    bestThreshold = thresholdRange[bestIdx]

    lowThreshold, highThreshold = chi2_confidence_interval(thresholdRange, chi2Curve, bestIdx)

    bestFitRatios = counts_to_ratios(get_multiplicity_counts(bestThreshold)[0])
    bestFitRateFull = seitz_count(bestFitRatios, sum(g["backSubFull"]))

    # same scan under the other normalization mode, purely for the side-by-side
    # comparison plot below - does not affect the best-fit threshold used elsewhere
    altChi2Curve = [chi_squared_calc(g, threshold, useGlobalNorm=not useGlobalChi2Normalization)
                     for threshold in thresholdRange]
    altBestIdx = int(np.argmin(altChi2Curve))
    altBestThreshold = thresholdRange[altBestIdx]
    altLowThreshold, altHighThreshold = chi2_confidence_interval(thresholdRange, altChi2Curve, altBestIdx)

    g["bestFit"] = {
        "thresholdRange": thresholdRange,
        "chi2Curve": chi2Curve,
        "threshold": bestThreshold,
        "thresholdErrLow": bestThreshold - lowThreshold,
        "thresholdErrHigh": highThreshold - bestThreshold,
        "chi2": chi2Curve[bestIdx],
        "rate": rebin(bestFitRateFull),
        "rateFull": bestFitRateFull,
        "altChi2Curve": altChi2Curve,
        "altThreshold": altBestThreshold,
        "altThresholdErrLow": altBestThreshold - altLowThreshold,
        "altThresholdErrHigh": altHighThreshold - altBestThreshold,
        "altChi2": altChi2Curve[altBestIdx],
    }

# chi-squared vs threshold scan
def plot_chi2_scan(gridKev, chi2Curve, bestThreshold, bestChi2, savepath):
    plt.figure(figsize=(8, 6))
    plt.plot(gridKev, chi2Curve, 'o', markersize=3, color="steelblue")
    plt.axvline(bestThreshold, color='red', linestyle='--',
                label=f"Best fit: {bestThreshold:0.2f} keV (" + r"$\chi^2$" + f"={bestChi2:0.2f})")
    plt.xlabel("Threshold [keV]", fontsize=16)
    plt.ylabel(r"$\chi^2$", fontsize=16)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig(savepath)
    plt.close()

# same scan, overlaying both normalization modes: blue is the mode currently selected
# by useGlobalChi2Normalization, red is the other one
def plot_chi2_scan_comparison(gridKev, chi2Curve, bestThreshold, bestChi2,
                               altChi2Curve, altBestThreshold, altBestChi2, savepath):
    label = normalization_mode_label(useGlobalChi2Normalization)
    altLabel = normalization_mode_label(not useGlobalChi2Normalization)

    plt.figure(figsize=(8, 6))
    plt.plot(gridKev, chi2Curve, 'o', markersize=3, color="steelblue", label=label)
    plt.axvline(bestThreshold, color='steelblue', linestyle='--',
                label=f"{label} best fit: {bestThreshold:0.2f} keV (" + r"$\chi^2$" + f"={bestChi2:0.2f})")
    plt.plot(gridKev, altChi2Curve, 'o', markersize=3, color="red", label=altLabel)
    plt.axvline(altBestThreshold, color='red', linestyle='--',
                label=f"{altLabel} best fit: {altBestThreshold:0.2f} keV (" + r"$\chi^2$" + f"={altBestChi2:0.2f})")
    plt.xlabel("Threshold [keV]", fontsize=16)
    plt.ylabel(r"$\chi^2$", fontsize=16)
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.savefig(savepath)
    plt.close()

for g in groups:
    plot_chi2_scan(
        g["bestFit"]["thresholdRange"], g["bestFit"]["chi2Curve"],
        g["bestFit"]["threshold"], g["bestFit"]["chi2"],
        savepath=plot_path(f"chi2Scans/chi2scan{g['p']}{g['T']}.png"),
    )
    plot_chi2_scan_comparison(
        g["bestFit"]["thresholdRange"], g["bestFit"]["chi2Curve"],
        g["bestFit"]["threshold"], g["bestFit"]["chi2"],
        g["bestFit"]["altChi2Curve"], g["bestFit"]["altThreshold"], g["bestFit"]["altChi2"],
        savepath=plot_path(f"chi2Scans/chi2scanCompare{g['p']}{g['T']}.png"),
    )

def plot_fit_to_seitz_ratio(groups, savepath):
    seitzVals = [g["seitz"] for g in groups]
    fitVals = [g["bestFit"]["threshold"] for g in groups]
    fitErrLow = [g["bestFit"]["thresholdErrLow"] for g in groups]
    fitErrHigh = [g["bestFit"]["thresholdErrHigh"] for g in groups]
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

plot_fit_to_seitz_ratio(groups, savepath=plot_path("fitToSeitzRatio.png"))

# if true, show delta-chi2=1 error bars on the comparison plot below
showComparisonErrorBars = True

# same plot, overlaying both normalization modes: blue is the mode currently selected
# by useGlobalChi2Normalization, red is the other one
def plot_fit_to_seitz_ratio_comparison(groups, savepath):
    label = normalization_mode_label(useGlobalChi2Normalization)
    altLabel = normalization_mode_label(not useGlobalChi2Normalization)

    seitzVals = [g["seitz"] for g in groups]
    fitVals = [g["bestFit"]["threshold"] for g in groups]
    altFitVals = [g["bestFit"]["altThreshold"] for g in groups]

    plt.figure(figsize=(8, 6))
    if showComparisonErrorBars:
        fitErrLow = [g["bestFit"]["thresholdErrLow"] for g in groups]
        fitErrHigh = [g["bestFit"]["thresholdErrHigh"] for g in groups]
        altFitErrLow = [g["bestFit"]["altThresholdErrLow"] for g in groups]
        altFitErrHigh = [g["bestFit"]["altThresholdErrHigh"] for g in groups]
        plt.errorbar(seitzVals, fitVals, yerr=[fitErrLow, fitErrHigh],
                     fmt='o', color="steelblue", ecolor="steelblue", capsize=4, label=label)
        plt.errorbar(seitzVals, altFitVals, yerr=[altFitErrLow, altFitErrHigh],
                     fmt='o', color="red", ecolor="red", capsize=4, label=altLabel)
    else:
        plt.plot(seitzVals, fitVals, 'o', color="steelblue", label=label)
        plt.plot(seitzVals, altFitVals, 'o', color="red", label=altLabel)
    plt.axline((0, 0), slope=1, color='gray', linestyle='--', label="Fit = Seitz")
    plt.xlabel("Seitz Threshold [keV]", fontsize=16)
    plt.ylabel("Best Fit Threshold [keV]", fontsize=16)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig(savepath)
    plt.close()

plot_fit_to_seitz_ratio_comparison(groups, savepath=plot_path("fitToSeitzRatioCompare.png"))

"""
SINLGE THRESHOLD RATES WITH AND WITHOUT THEORETICAL THRESHOLD
"""
binLabels3 = ["1", "2", "3+"]
binLabelsFull = ["1", "2", "3", "4", "5+"]

zeroKevCountsRaw, _ = get_multiplicity_counts(0.0)
zeroKevCountsRebinned = rebin(zeroKevCountsRaw)
zeroKevRate3 = [normalizationFactor * c / sum(zeroKevCountsRebinned) for c in zeroKevCountsRebinned]
zeroKevRatiosFull = counts_to_ratios(zeroKevCountsRaw)

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

for g in groups:
    if useRebinnedThresholdPlots:
        linhistArgs = (
            binLabels3, rebin(g["binCounts"]), rebin_errors(g["binCountError"]),
            rebin(g["backBins"]), rebin_errors(g["backErrorLow"]), rebin_errors(g["backErrorHigh"]),
            g["backSub"], g["errLow"], g["errHigh"],
            zeroKevRate3, g["seitzRate"], g["seitz"],
        )
        bestFitRate = g["bestFit"]["rate"]
    else:
        total = sum(g["backSubFull"])
        linhistArgs = (
            binLabelsFull, g["binCounts"], g["binCountError"], g["backBins"], g["backErrorLow"], g["backErrorHigh"],
            g["backSubFull"], g["backSubErrorLowFull"], g["backSubErrorHighFull"],
            seitz_count(zeroKevRatiosFull, total), seitz_count(counts_to_ratios(g["seitzCountsFull"]), total), g["seitz"],
        )
        bestFitRate = g["bestFit"]["rateFull"]

    # one version without the chi-squared best-fit curve, one with it
    plot_linhist(*linhistArgs, savepath=plot_path(f"linHists/linhist{g['p']}{g['T']}.png"))
    plot_linhist(
        *linhistArgs, savepath=plot_path(f"linHists/linhistFit{g['p']}{g['T']}.png"),
        bestFitRate=bestFitRate, bestThreshold=g["bestFit"]["threshold"], bestChi2=g["bestFit"]["chi2"],
    )


"""
AVERAGE THRESHOLD RATE WITH AND WITHOUT THEORETICAL THRESHOLD
"""
# bascially just redo everything but with all the data
binCountsAvg, sourceTimeAvg = bin_multiplicities(
    bubbleCount, sourceTimes, keep=lambda i, region: region not in excludedRegions
)
backBinsAvg, backErrorLowAvg, backErrorHighAvg, backSubBinsAvg, backSubErrorLowAvg, backSubErrorHighAvg, binCountErrorAvg = background_subtract(
    binCountsAvg, sourceTimeAvg, backgroundBinCounts, backgroundTime
)

sourceTimeAvg /= 60
binCountsAvg = [v / sourceTimeAvg for v in binCountsAvg]
binCountErrorAvg = [v / sourceTimeAvg for v in binCountErrorAvg]
backBinsAvg = [v / sourceTimeAvg for v in backBinsAvg]
backErrorLowAvg = [v / sourceTimeAvg for v in backErrorLowAvg]
backErrorHighAvg = [v / sourceTimeAvg for v in backErrorHighAvg]
backSubBinsAvg = [v / sourceTimeAvg for v in backSubBinsAvg]
backSubErrorLowAvg = [v / sourceTimeAvg for v in backSubErrorLowAvg]
backSubErrorHighAvg = [v / sourceTimeAvg for v in backSubErrorHighAvg]

avgSeitz = np.mean([g["seitz"] for g in groups])
avgSeitzCountsRaw, _ = get_multiplicity_counts(avgSeitz)

avgGroup = {
    "liveTime": sourceTimeAvg,
    "backSub": rebin(backSubBinsAvg),
    "errLow": rebin_errors(backSubErrorLowAvg),
    "errHigh": rebin_errors(backSubErrorHighAvg),
    "backSubFull": backSubBinsAvg,
    "backSubErrorLowFull": backSubErrorLowAvg,
    "backSubErrorHighFull": backSubErrorHighAvg,
}
avgChi2Curve = [chi_squared_calc(avgGroup, threshold) for threshold in thresholdRange]
avgBestIdx = int(np.argmin(avgChi2Curve))
avgBestThreshold = thresholdRange[avgBestIdx]
avgBestFitRatios = counts_to_ratios(get_multiplicity_counts(avgBestThreshold)[0])
avgBestFitRateFull = seitz_count(avgBestFitRatios, sum(backSubBinsAvg))

avgAltChi2Curve = [chi_squared_calc(avgGroup, threshold, useGlobalNorm=not useGlobalChi2Normalization)
                    for threshold in thresholdRange]
avgAltBestIdx = int(np.argmin(avgAltChi2Curve))
avgAltBestThreshold = thresholdRange[avgAltBestIdx]

avgBestFit = {
    "thresholdRange": thresholdRange,
    "chi2Curve": avgChi2Curve,
    "threshold": avgBestThreshold,
    "chi2": avgChi2Curve[avgBestIdx],
    "rate": rebin(avgBestFitRateFull),
    "rateFull": avgBestFitRateFull,
    "altChi2Curve": avgAltChi2Curve,
    "altThreshold": avgAltBestThreshold,
    "altChi2": avgAltChi2Curve[avgAltBestIdx],
}
plot_chi2_scan(
    avgBestFit["thresholdRange"], avgBestFit["chi2Curve"],
    avgBestFit["threshold"], avgBestFit["chi2"],
    savepath=plot_path("chi2scanAvg.png"),
)
plot_chi2_scan_comparison(
    avgBestFit["thresholdRange"], avgBestFit["chi2Curve"],
    avgBestFit["threshold"], avgBestFit["chi2"],
    avgBestFit["altChi2Curve"], avgBestFit["altThreshold"], avgBestFit["altChi2"],
    savepath=plot_path("chi2scanCompareAvg.png"),
)

if useRebinnedThresholdPlots:
    avgSeitzCountsRebinned = rebin(avgSeitzCountsRaw)
    avgSeitzRate = [normalizationFactor * c / sum(avgSeitzCountsRebinned) for c in avgSeitzCountsRebinned]
    avgLinhistArgs = (
        binLabels3, rebin(binCountsAvg), rebin_errors(binCountErrorAvg),
        rebin(backBinsAvg), rebin_errors(backErrorLowAvg), rebin_errors(backErrorHighAvg),
        rebin(backSubBinsAvg), rebin_errors(backSubErrorLowAvg), rebin_errors(backSubErrorHighAvg),
        zeroKevRate3, avgSeitzRate, avgSeitz,
    )
    avgBestFitRate = avgBestFit["rate"]
else:
    totalAvg = sum(backSubBinsAvg)
    avgSeitzRatiosFull = counts_to_ratios(avgSeitzCountsRaw)
    avgLinhistArgs = (
        binLabelsFull, binCountsAvg, binCountErrorAvg, backBinsAvg, backErrorLowAvg, backErrorHighAvg,
        backSubBinsAvg, backSubErrorLowAvg, backSubErrorHighAvg,
        seitz_count(zeroKevRatiosFull, totalAvg), seitz_count(avgSeitzRatiosFull, totalAvg), avgSeitz,
    )
    avgBestFitRate = avgBestFit["rateFull"]

# one version without the chi-squared best-fit curve, one with it
plot_linhist(*avgLinhistArgs, savepath=plot_path("avgseitz.png"))
plot_linhist(
    *avgLinhistArgs, savepath=plot_path("avgseitzFit.png"),
    bestFitRate=avgBestFitRate, bestThreshold=avgBestFit["threshold"], bestChi2=avgBestFit["chi2"],
)


