import matplotlib.pyplot as plt
import numpy as np
import sys, csv, glob, os
from sbcbinaryformat import Streamer, Writer

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
backgroundList = []
for run in backgroundRunsWarm:
    backgroundList.append(run)
for run in backgroundRunsCold:
    backgroundList.append(run)
for run in backgroundRunsHot:
    backgroundList.append(run)

def process_dir_background(dirpath):
    checkedRuns = []
    returnList = []
    for path in glob.glob(os.path.join(dirpath, "*.txt")):
        with open(path, 'r', encoding='utf-8') as f:
            for raw in f:
                parts = raw.split()
                if len(parts) < 5:
                    continue
                run = parts[0]
                ev = parts[1]
                try:
                    mult = int(parts[4])
                    region =int(parts[3])
                except ValueError:
                    continue
                if run not in backgroundList:
                    checkedRuns.append((ev,run))
                    continue
                if (ev,run) in checkedRuns:
                    continue
                returnList.append((run,ev,mult))
                checkedRuns.append((ev,run))
    return returnList

handScannedBackgrounds = process_dir_background('/exp/e961/data/SBC-25-handscan/')



def process_dir_ana(dirpath):
    outList = []

    subdirs = {os.path.basename(p.rstrip(os.sep)) for p in glob.glob(os.path.join(dirpath, '*/'))}
    for item in handScannedBackgrounds:
        dirname = item[0]
        if dirname in subdirs:
            fullPath = os.path.join(dirpath, dirname)
            expData = Streamer(fullPath + '/exposure.sbc').to_dict()
            quickCheck = False
            for i in range(len(expData["ev"])):
                if int(expData["ev"][i]) == int(item[1]) and float(expData['PT2121_livetime'][i]) > float(1):
                    outList.append((int(item[2]), float(expData['PT2121_livetime'][i])))
                    quickCheck = True
                    break

    return outList



backgroundPairs = process_dir_ana('/exp/e961/data/SBC-25-recon/dev-output/')
backgroundTime = 0
backgroundSingles = 0
background2s = 0
background3s = 0
background4s = 0
background5s = 0

for i in range(len(backgroundPairs)):
    if backgroundPairs[i][0] == 1:
        backgroundSingles += 1
    elif backgroundPairs[i][0] == 2:
        background2s +=1
    elif backgroundPairs[i][0] == 3:
        background3s +=1
    elif backgroundPairs[i][0] == 4:
        background4s +=1
    elif backgroundPairs[i][0] >= 5:
        background5s +=1
    backgroundTime += backgroundPairs[i][1]
backgroundMultis = background2s + background3s + background4s + background5s

# neutron stuff

# path to csv output from handscanning in LAr10Ana
path = "/exp/e961/data/SBC-25-handscan/"

## warm annular
neutronRunsWarm = ["20260107_3", "20260107_4", "20260107_5", "20260107_6", "20260107_7", "20260108_0", "20260108_1", "20260108_2", "20260108_3", "20260108_4", "20260108_5", "20260108_6", "20260108_7", "20260108_8", "20260109_0"]
## cold annular
neutronRunsCold = ["20260122_3","20260122_4","20260122_5","20260122_6","20260123_0","20260123_1","20260123_2","20260123_3","20260123_4","20260123_8","20260123_9","20260123_10","20260124_0","20260124_1","20260124_3","20260124_4","20260124_5","20260125_0","20260125_1","20260125_2","20260125_3","20260125_4","20260125_5","20260125_6","20260125_7","20260125_8"]

## 119K
neutronRunsHot = ["20260205_12","20260205_13","20260205_14","20260205_15","20260205_16","20260205_17","20260205_18","20260206_0","20260206_1","20260206_2","20260206_3","20260206_4","20260206_5","20260206_6","20260206_7","20260213_1","20260213_2","20260213_3","20260213_4","20260213_5","20260213_6","20260213_7","20260213_8","20260213_9","20260214_0","20260214_1","20260214_2","20260214_3","20260214_4","20260214_5","20260214_6","20260214_7","20260214_8","20260214_9","20260214_10","20260214_11","20260214_12","20260214_13","20260214_14","20260215_0","20260215_1","20260215_2","20260215_3","20260215_4","20260215_5","20260215_6","20260215_7","20260215_8","20260215_9","20260215_10","20260215_11","20260215_12","20260215_13","20260215_14","20260216_0","20260216_1","20260216_2","20260216_3","20260216_4","20260216_5","20260216_6","20260216_7","20260216_8","20260216_9","20260216_10","20260216_11","20260216_12","20260216_13","20260217_0","20260217_1","20260217_2","20260217_3","20260217_4","20260217_5","20260217_6"]


## ones that are used for this graph
neutronRuns = []
for run in neutronRunsWarm:
    neutronRuns.append(run)
for run in neutronRunsCold:
    neutronRuns.append(run)
for run in neutronRunsHot:
    neutronRuns.append(run)

#neutronRuns.append("20260330_5")
## purley for the sake of me testing

def process_dir(dirpath):
    checkedRuns = []
    bubbleCount = []
    sourceTime = 0
    for path in glob.glob(os.path.join(dirpath, "*.txt")):
        with open(path, 'r', encoding='utf-8') as f:
            for raw in f:
                parts = raw.split()
                if len(parts) < 5:
                    continue
                run = parts[0]
                ev = parts[1]
                try:
                    mult = int(parts[4])
                    region =int(parts[3])
                except ValueError:
                    continue
                if run not in neutronRuns:
                    checkedRuns.append((ev,run))
                    continue
                if (ev,run) in checkedRuns:
                    continue
                bubbleCount.append((mult,region))
                checkedRuns.append((ev,run))
                expData = Streamer('/exp/e961/data/SBC-25-recon/dev-output/' + run + '/exposure.sbc').to_dict()
                for i in range(len(expData["ev"])):
                    if int(expData["ev"][i]) == int(ev):
                        sourceTime += float(expData['PT2121_livetime'][i])
                        break
    return bubbleCount, sourceTime

bubbleCount, sourceTime = process_dir(path)

# bin data and calc ratios 
binLabels = ["1", "2", "3", "4", "5+"]
binCounts = [0]*(len(binLabels))
# idk if we need this but could help
excludedRegions = []
for n in bubbleCount:
    if n[1] in excludedRegions:
        continue
    if n[0] >= 5:
            binCounts[4] += 1
    else:
        binCounts[n[0]-1] += 1
# background rate subtraction
backgroundSingleEst = backgroundSingles * sourceTime/backgroundTime 
background2sEst = background2s * sourceTime/backgroundTime 
background3sEst = background3s * sourceTime/backgroundTime 
background4sEst = background4s * sourceTime/backgroundTime 
background5sEst = background5s * sourceTime/backgroundTime 
backBins = [backgroundSingleEst, background2sEst, background3sEst, background4sEst, background5sEst]
backSubBins = []
for i in range(len(binCounts)):
    backSubBins.append(binCounts[i] - backBins[i])


binCountError = []
binCountError.append(np.sqrt(0))
ratios = [1]
backSubRatios = [1]
ratioError = []
for c in binCounts[1:]:
    ratios.append(c/binCounts[0])
    binCountError.append(np.sqrt(c))    
for c in backSubBins[1:]:
    backSubRatios.append(c/backSubBins[0])


for i in range(0,len(binLabels)):
    ratioError.append( np.sqrt(np.abs( (binCountError[i]/binCounts[0])**2 + (binCounts[i] * binCountError[0]/(binCounts[0]**2))**2   )))


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

simCountMin = []
simCountMax = []
for i in range(0,len(thresholds)-1):
    simCountMin.append([])
    simCountMax.append([])
    for j in range(0,len(ratiosSim[i])-1):
        simCountMin[i].append(np.abs(binCounts[0]* (ratiosSim[j][i] - simError[j][i])))
        simCountMax[i].append(np.abs(binCounts[0] * (ratiosSim[j][i] + simError[j][i])))






# graphing
x = np.arange(len(binLabels))

plt.figure(figsize=(16,9))
barsList = []
width = 0.9/(len(thresholds)-1)
colors = ["blue","red","green","orange", "teal","black"]
num_groups = len(thresholds) - 1

for i in range(num_groups):
    offset = (i - (num_groups - 1) / 2) * width
    xPos = x + offset
    bars = plt.bar(xPos, simCountMax[i], width=width, color=colors[i % len(colors)],
                   edgecolor="black", alpha=0.18,zorder=0)
    barsList.append(bars)


for i in range(num_groups):
    for j, bar in enumerate(barsList[i]):
        r = ratiosSim[j][i]    # note swapped indices
        cx = bar.get_x() + bar.get_width()/2
        cy = bar.get_height()
        plt.text(cx, cy + 0.01*max(binCounts), f"{r:0.3f}", ha='center', va='bottom', fontsize=12)

shadedInBars = []
for i in range(num_groups):
    offset = (i - (num_groups - 1) / 2) * width
    xPos = x + offset
    bars = plt.bar(xPos, simCountMin[i], width=width, color=colors[i % len(colors)],
                   edgecolor="black", label=f"{thresholds[i]/1000}KeV",zorder=2)
    shadedInBars.append(bars)

points = x
plt.errorbar(points, binCounts, yerr=binCountError,fmt='o',color="red", ecolor="red", label="Source Rate")

# subtracted rates
plt.errorbar(points, backBins,fmt='o',color="orange", label="Background Rate")
plt.errorbar(points, backSubBins,fmt='o',color="purple", label="Background Subtracted Rate")


plt.xticks(x,binLabels)
plt.xlabel("Bubble Multiplicity",fontsize=16)
plt.ylabel("Count",fontsize=16)
plt.title("Handscan Bubble Multiplicites and Simulation Comparison",fontsize=16)
plt.legend(title="Thresholds",fontsize=16)
plt.savefig("multhist.png")
plt.show()



