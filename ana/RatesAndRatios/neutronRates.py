import matplotlib.pyplot as plt
import numpy as np
import glob, os
from sbcbinaryformat import Streamer
import SeitzModel as sm

# background rate caluclation for subtraction
backgroundRunsWarm = [
"20251113_9",
"20251113_10",
"20251113_11",
"20251114_0",
"20251114_1",
"20251114_6",
"20251114_36",
"20251114_37",
"20251115_0",
"20251115_1",
"20251115_2",
"20251115_3",
"20251115_4",
"20251115_5",
"20251116_1",
"20251116_2",
"20251117_0",
"20251117_1",
"20251126_7",
"20251126_8",
"20251127_0",
"20251127_1",
"20251127_2",
"20251127_3",
"20251127_4",
"20251127_5",
"20251128_0",
"20251128_1",
"20251128_2",
"20251128_3",
"20251128_4",
"20251129_0",
"20251129_1",
"20251129_2",
"20251129_3",
"20251129_4",
"20251129_5",
"20251130_0",
"20251130_1",
"20251130_2",
"20251130_3",
"20251130_4",
"20251130_5",
]
## cold annular
backgroundRunsCold = [
"20260117_0",
"20260117_1",
"20260117_2",
"20260117_3",
"20260117_4",
"20260118_0",
"20260118_1",
"20260118_2",
"20260118_3",
"20260118_4",
"20260119_0",
"20260119_1",
"20260119_2",
"20260120_0",
"20260120_1",
]
## 199k
backgroundRunsHot = [
"20260217_7",
"20260217_8",
"20260217_9",
"20260217_10",
"20260217_11",
"20260217_12",
"20260217_13",
"20260218_0",
"20260218_1",
"20260218_2",
"20260218_3",
"20260218_4",
"20260218_5",
"20260218_6",
"20260218_15",
"20260218_16",
"20260219_0",
"20260219_1",
"20260219_2",
"20260219_3",
"20260219_4",
"20260219_5",
"20260219_6",
"20260219_7",
"20260219_8",
"20260219_9",
"20260219_10",
"20260219_11",
"20260220_1",
"20260220_2",
"20260220_3",
"20260220_4",
]


## ones to use for rate calculation
backgroundList = backgroundRunsWarm + backgroundRunsCold + backgroundRunsHot


# yields (run, ev, mult, region) once per unique event whose run is in runList
# shared by process_dir_background and process_dir since they were doing the same parsing
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


def process_dir_background(dirpath):
    return [(run, ev, mult) for run, ev, mult, _ in iter_matched_events(dirpath, backgroundList)]

handScannedBackgrounds = process_dir_background('/exp/e961/data/SBC-25-handscan/')


def process_dir_ana(dirpath):
    outList = []
    subdirs = {os.path.basename(p.rstrip(os.sep)) for p in glob.glob(os.path.join(dirpath, '*/'))}
    for run, ev, mult in handScannedBackgrounds:
        if run not in subdirs:
            continue
        expData = Streamer(os.path.join(dirpath, run, 'exposure.sbc')).to_dict()
        for i in range(len(expData["ev"])):
            if int(expData["ev"][i]) == int(ev) and float(expData['PT2121_livetime'][i]) > 1:
                outList.append((mult, float(expData['PT2121_livetime'][i])))
                break
    return outList


backgroundPairs = process_dir_ana('/exp/e961/data/SBC-25-recon/dev-output/')

backgroundBinCounts = [0] * 5
backgroundTime = 0.0
for mult, livetime in backgroundPairs:
    backgroundBinCounts[min(mult, 5) - 1] += 1
    backgroundTime += livetime
backgroundSingles, background2s, background3s, background4s, background5s = backgroundBinCounts
backgroundMultis = background2s + background3s + background4s + background5s

# neutron stuff

# path to csv output from handscanning in LAr10Ana
path = "/exp/e961/data/SBC-25-handscan/"



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

## the multhist bar plots are slow and not always wanted, toggle off when not needed
makeMultPlots = True

## whether the background-subtracted rate's error bar carries through the background
## error's asymmetry (zero-count bins only have an upper limit), or stays symmetric
asymmetricBackSubError = False


def process_dir(dirpath):
    bubbleCount = []
    sourceTimes = []
    ptemps = []
    for run, ev, mult, region in iter_matched_events(dirpath, neutronRuns):
        bubbleCount.append((mult, region))
        expData = Streamer(f'/exp/e961/data/SBC-25-recon/dev-output/{run}/exposure.sbc').to_dict()
        for i in range(len(expData["ev"])):
            if int(expData["ev"][i]) == int(ev):
                sourceTimes.append(float(expData['PT2121_livetime'][i]))
                break
        evData = Streamer(f'/exp/e961/data/SBC-25-recon/dev-output/{run}/event.sbc').to_dict()
        for i in range(len(evData["ev"])):
            if int(evData["ev"][i]) == int(ev):
                # NOTE: preserved as-is; this checks neutronRunsHot (config A, always empty)
                # rather than neutronRunsHotB, so T is always 116 here. Flagging, not fixing.
                temp = 119 if run in neutronRunsHot else 116
                ptemps.append((evData["pset_lo"][i], evData["pset_hi"][i], temp))
                break
    return bubbleCount, sourceTimes, ptemps

bubbleCount, sourceTimes, psetsTemps  = process_dir(path)
pToUse = []
for psets in psetsTemps:
    if float(psets[0]) == float(psets[1]) and (float(psets[0]), float(psets[2])) not in pToUse:
        pToUse.append((float(psets[0]), float(psets[2])))
pToUse.sort()

# this is probably not the best way to do this but it works for now
# thresholds in eV, case B from ryan
thresholds = [0.0, 5000.0, 10000.0, 15000.0, 20000.0, 25000.0]
ratiosSim = []
simError  = []
ratiosSim.append([1,1,1,1,1,1])
simError.append([0,0,0,0,0,0])
ratiosSim.append([0.40904175, 0.3605948, 0.33270232,  0.31065533, 0.2928382,  0.26341731])
simError.append([0.01542234, 0.0195329, 0.0197859, 0.01941437, 0.01920887, 0.01817524])
ratiosSim.append([0.18451222, 0.12556795, 0.10411737, 0.08454227, 0.07374005, 0.06626506])
simError.append([0.00903129, 0.00975382, 0.00928462, 0.00839415, 0.00795298, 0.00757482])
ratiosSim.append([0.08392819, 0.04750103, 0.03265499, 0.02751376, 0.02068966, 0.01642935])
simError.append([0.00549518, 0.00539489, 0.00464159, 0.00432533, 0.00378953, 0.00338405])
ratiosSim.append([0.03417694, 0.01363073, 0.01183152, 0.00800400, 0.00636605, 0.00657174])
simError.append([0.00322163, 0.00264983, 0.00262369, 0.00218002, 0.00198435, 0.00205089])

seitzRatios = []
seitzThresholds = [1.0530287933286058, 1.1391899070606315, 1.235873456434369, 1.3448015203144668, 1.468053352541133, 1.7682026967519062, 2.164316091009131, 2.699825488759606, 3.4446953964301947, 4.516768136484239, 8.668378138755212, 20.90193280201993]

seitzRatios.append([1,1,1,1,1,1,1,1,1,1,1,1])
seitzRatios.append([0.4151436,0.41253264,0.41469816,0.41315789,0.41052632,0.40789474,0.40419948,0.4047619,0.40909091,0.41734417,0.39678284,0.38419619])
seitzRatios.append([0.14882507,0.14882507,0.1496063,0.15,0.15,0.15,0.1496063,0.15079365,0.15240642,0.14905149,0.14477212,0.14945652])
seitzRatios.append([0.08616188,0.08616188,0.08661417,0.08684211,0.08684211,0.08684211,0.08661417,0.08994709,0.09090909,0.09214092,0.09115282,0.08991826])
seitzRatios.append([0.02610966,0.02610966,0.02624672,0.02631579,0.02631579,0.02631579,0.02624672,0.02380952,0.02406417,0.02439024,0.02412869,0.02724796])

# bin data and calc ratios
# keep(i, region) picks which events count towards this histogram
def bin_multiplicities(bubbleCount, sourceTimes, keep):
    binCounts = [0] * 5
    sourceTime = 0.0
    count = 0
    for i, (mult, region) in enumerate(bubbleCount):
        if not keep(i, region):
            continue
        sourceTime += sourceTimes[i]
        count += 1
        binCounts[min(mult, 5) - 1] += 1
    return binCounts, sourceTime, count


# background rate subtraction
# make sure this is in a rate
def background_subtract(binCounts, sourceTime, backgroundBinCounts, backgroundTime, asymmetricBackSubError):
    backBins = [b * sourceTime / backgroundTime for b in backgroundBinCounts]
    # b < 1 -> sqrt(b) alone would push the lower bound negative, so switch to a
    # one-sided upper limit there too, not just at exactly 0. this has to be checked
    # on backBins (before the caller scales everything to a per-minute rate) - if
    # checked after, every bin's rate is small and this would trigger everywhere
    backErrorLow = [np.sqrt(b) if b >= 1 else 0.0 for b in backBins]
    backErrorHigh = [np.sqrt(b) if b >= 1 else -np.log(1 - 0.68) / backgroundTime for b in backBins]
    binCountError = [np.sqrt(c) for c in binCounts]
    backSubBins = [c - b for c, b in zip(binCounts, backBins)]
    # background high -> subtracted rate pulled down; background low -> subtracted rate pulled up
    backSubErrorLow = [np.sqrt(np.abs(ce**2 + be**2)) for ce, be in zip(binCountError, backErrorHigh)]
    if asymmetricBackSubError:
        backSubErrorHigh = [np.sqrt(np.abs(ce**2 + be**2)) for ce, be in zip(binCountError, backErrorLow)]
    else:
        # must be an independent list, not an alias: the caller divides both lists
        # element-wise by sourceTime later, which would double-divide a shared list
        backSubErrorHigh = list(backSubErrorLow)
    return backBins, backErrorLow, backErrorHigh, backSubBins, backSubErrorLow, backSubErrorHigh, binCountError


# this is probably not the best way to do this but it works for now
# normalize each threshold's simulated multiplicity shape so its total matches
# the background-subtracted total (total), instead of anchoring to a single bin
def sim_count_bounds(thresholds, ratiosSim, simError, total):
    simCountMin, simCountMax = [], []
    for i in range(len(thresholds) - 1):
        ratios_i = [ratiosSim[j][i] for j in range(len(ratiosSim))]
        errors_i = [simError[j][i] for j in range(len(ratiosSim))]
        scale = total / sum(ratios_i)
        simCountMin.append([np.abs(scale * (r - e)) for r, e in zip(ratios_i, errors_i)])
        simCountMax.append([np.abs(scale * (r + e)) for r, e in zip(ratios_i, errors_i)])
    return simCountMin, simCountMax


# seitz threshold
# same idea: scale the simulated shape at the seitz threshold to match the
# background-subtracted total, instead of anchoring to backSubBins[0]
def seitz_count(seitz, seitzThresholds, seitzRatios, thresholds, ratiosSim, total):
    if seitz in seitzThresholds:
        col = seitzThresholds.index(seitz)
        ratios = [seitzRatios[j][col] for j in range(len(seitzRatios))]
    else:
        print("Couldnt find exact value, using lerp")
        # linear interpolate for now
        ratios = None
        for i in range(1, len(thresholds)):
            if thresholds[i - 1] <= seitz <= thresholds[i]:
                t = (seitz - thresholds[i - 1]) / (thresholds[i] - thresholds[i - 1])
                ratios = [ratiosSim[j][i - 1] + (ratiosSim[j][i] - ratiosSim[j][i - 1]) * t for j in range(len(ratiosSim))]
                break
    scale = total / sum(ratios)
    return [scale * r for r in ratios]


def plot_multiplicity_bars(binLabels, thresholds, ratiosSim, simCountMin, simCountMax,
                            binCounts, binCountError, backBins, backErrorLow, backErrorHigh,
                            backSubBins, backSubErrorLow, backSubErrorHigh,
                            title, savepath):
    x = np.arange(len(binLabels))
    width = 0.9 / (len(thresholds) - 1)
    colors = ["blue", "red", "green", "orange", "teal", "black"]
    num_groups = len(thresholds) - 1

    plt.figure(figsize=(16, 9))
    barsList = []
    for i in range(num_groups):
        offset = (i - (num_groups - 1) / 2) * width
        bars = plt.bar(x + offset, simCountMax[i], width=width, color=colors[i % len(colors)],
                        edgecolor="black", alpha=0.18, zorder=0)
        barsList.append(bars)

    for i in range(num_groups):
        for j, bar in enumerate(barsList[i]):
            r = ratiosSim[j][i]    # note swapped indices
            cx = bar.get_x() + bar.get_width() / 2
            cy = bar.get_height()
            plt.text(cx, cy + 0.01 * max(binCounts), f"{r:0.3f}", ha='center', va='bottom', fontsize=12)

    for i in range(num_groups):
        offset = (i - (num_groups - 1) / 2) * width
        plt.bar(x + offset, simCountMin[i], width=width, color=colors[i % len(colors)],
                edgecolor="black", label=f"{thresholds[i]/1000}keV", zorder=2)

    plt.errorbar(x, binCounts, yerr=binCountError, fmt='o', color="red", ecolor="red", label="Source Rate")

    # subtracted rates
    plt.errorbar(x, backBins, yerr=[backErrorLow, backErrorHigh], fmt='o', color="orange", label="Background Rate")
    plt.errorbar(x, backSubBins, yerr=[backSubErrorLow, backSubErrorHigh], fmt='o', color="purple", label="Background Subtracted Rate")

    plt.xticks(x, binLabels)
    plt.xlabel("Bubble Multiplicity", fontsize=20)
    plt.ylabel("Count", fontsize=20)
    plt.title(title, fontsize=20)
    plt.legend(title="Thresholds", fontsize=20)
    plt.savefig(savepath)
    plt.close()


binLabels = ["1", "2", "3", "4", "5+"]
#idk if we need this but could help
excludedRegions = []

usedSeitz = []
mostCommon = None
mostCommonCount = 0
for p, T in pToUse:
    binCounts, sourceTime, curCount = bin_multiplicities(
        bubbleCount, sourceTimes,
        keep=lambda i, region: region not in excludedRegions and psetsTemps[i][0] == p and psetsTemps[i][2] == T,
    )
    if curCount > mostCommonCount:
        mostCommonCount = curCount
        mostCommon = (p, T)

    # background rate subtraction
    # make sure this is in a rate
    backBins, backErrorLow, backErrorHigh, backSubBins, backSubErrorLow, backSubErrorHigh, binCountError = background_subtract(
        binCounts, sourceTime, backgroundBinCounts, backgroundTime, asymmetricBackSubError
    )

    # change this into a rate plot instead of count plot
    # change to minutes?
    sourceTime /= 60
    for lst in (binCounts, binCountError, backBins, backErrorLow, backErrorHigh,
                backSubBins, backSubErrorLow, backSubErrorHigh):
        for i in range(len(lst)):
            lst[i] /= sourceTime

    # normalize the simulated/seitz predictions to the background-subtracted total
    meanNorm = sum(backSubBins)
    simCountMin, simCountMax = sim_count_bounds(thresholds, ratiosSim, simError, meanNorm)

    # seitz threshold
    seitz = sm.SeitzModel(p * 14.5038, -273.15 + T, 'argon').Q
    usedSeitz.append(seitz)
    seitzCount = seitz_count(seitz, seitzThresholds, seitzRatios, thresholds, ratiosSim, meanNorm)

    # new graphing
    plt.figure(figsize=(10, 10))
    x = np.arange(len(binLabels))
    edges = np.concatenate(([x[0] - 0.5], (x[:-1] + x[1:]) / 2, [x[-1] + 0.5]))

    # source on
    plt.errorbar(x, binCounts, yerr=binCountError, fmt='o', color="red", ecolor="red", label="Source Rate")
    # source off
    plt.errorbar(x, backBins, yerr=[backErrorLow, backErrorHigh], fmt='o', color="blue", label="Background Rate")
    # on - off
    plt.errorbar(x, backSubBins, yerr=[backSubErrorLow, backSubErrorHigh], fmt='o', color="purple", label="Background Subtracted Rate")

    # 0KeV threshold
    plt.stairs(simCountMax[0], edges, color="orange", linewidth=4, label="0keV Threshold")
    plt.stairs(seitzCount, edges, color="green", linewidth=6, label=f"Seitz Threshold\n({seitz:0.2f} keV )")

    plt.xticks(x, binLabels, fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlabel("Bubble Multiplicity", fontsize=20)
    plt.ylabel("Rate [count/min]", fontsize=20)
    plt.legend(fontsize=20)
    plt.tight_layout()
    plt.savefig(f"linhist{p}{T}.png")
    plt.close()

    # old graphing
    if makeMultPlots:
        plot_multiplicity_bars(
            binLabels, thresholds, ratiosSim, simCountMin, simCountMax,
            binCounts, binCountError, backBins, backErrorLow, backErrorHigh,
            backSubBins, backSubErrorLow, backSubErrorHigh,
            title=f"Multiplicites Comparison for P={p}bara and T={T}K",
            savepath=f"multhist{p}{T}.png",
        )

# averaged out stats
binCounts, sourceTime, _ = bin_multiplicities(
    bubbleCount, sourceTimes, keep=lambda i, region: region not in excludedRegions
)
backBins, backErrorLow, backErrorHigh, backSubBins, backSubErrorLow, backSubErrorHigh, binCountError = background_subtract(
    binCounts, sourceTime, backgroundBinCounts, backgroundTime, asymmetricBackSubError
)

# change this into a rate plot instead of count plot
# change to minutes?
sourceTime /= 60
for lst in (binCounts, binCountError, backBins, backErrorLow, backErrorHigh,
            backSubBins, backSubErrorLow, backSubErrorHigh):
    for i in range(len(lst)):
        lst[i] /= sourceTime

# normalize the simulated predictions to the background-subtracted total
meanNorm = sum(backSubBins)
simCountMin, simCountMax = sim_count_bounds(thresholds, ratiosSim, simError, meanNorm)

if makeMultPlots:
    plot_multiplicity_bars(
        binLabels, thresholds, ratiosSim, simCountMin, simCountMax,
        binCounts, binCountError, backBins, backErrorLow, backErrorHigh,
        backSubBins, backSubErrorLow, backSubErrorHigh,
        title="Multiplicites Comparison for all P and T",
        savepath="multhistavg.png",
    )

# averaged seitz threshold
avg = np.mean(usedSeitz)
if avg != 4.200422148324119:
    print("go get the new avg threshold it's value is:" + str(avg))
    exit()
ratios = [1, 0.41018767, 0.15013405, 0.09115282, 0.02412869]

plt.figure(figsize=(10, 10))
x = np.arange(len(binLabels))
edges = np.concatenate(([x[0] - 0.5], (x[:-1] + x[1:]) / 2, [x[-1] + 0.5]))

# source on
plt.errorbar(x, binCounts, yerr=binCountError, fmt='o', color="red", ecolor="red", label="Source Rate")
# source off
plt.errorbar(x, backBins, yerr=[backErrorLow, backErrorHigh], fmt='o', color="blue", label="Background Rate")
# on - off
plt.errorbar(x, backSubBins, yerr=[backSubErrorLow, backSubErrorHigh], fmt='o', color="purple", label="Background Subtracted Rate")

# 0KeV threshold
plt.stairs(simCountMax[0], edges, color="orange", linewidth=4, label="0keV Threshold")
# seitz threshold
# scale the averaged seitz shape to match the background-subtracted total, same as the per (P,T) loop
seitz = avg
seitzCount = [(meanNorm / sum(ratios)) * r for r in ratios]

plt.stairs(seitzCount, edges, color="green", linewidth=6, label=f"Seitz Threshold ({seitz:0.2f} keV)")

plt.xticks(x, binLabels, fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel("Bubble Multiplicity", fontsize=20)
plt.ylabel("Rate [count/min]", fontsize=20)
plt.legend(fontsize=20)
plt.tight_layout()
plt.savefig("avgseitz.png", dpi=500)
plt.close()


