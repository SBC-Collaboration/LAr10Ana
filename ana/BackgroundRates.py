import matplotlib.pyplot as plt
import numpy as np
import sys, csv, glob, os
from sbcbinaryformat import Streamer, Writer

## warm annular
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



def process_dir_txt(dirpath):
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


handScannedBackgrounds = process_dir_txt('/exp/e961/data/SBC-25-handscan/')

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

            if not quickCheck:
                print("never added")
    return outList

backgroundPairs = process_dir_ana('/exp/e961/data/SBC-25-recon/dev-output/')
print(backgroundPairs)
backgroundTime = 0
backgroundSingles = 0
backgroundMultis  = 0
for i in range(len(backgroundPairs)):
    if backgroundPairs[i][0] == 1:
        backgroundSingles += 1
    elif backgroundPairs[i][0] != 0:
        backgroundMultis += 1
    backgroundTime += backgroundPairs[i][1]

#convert to per hour
backgroundTime *= 1/(60*60)
print("Background sinlges rate is "+ str(backgroundSingles/backgroundTime) +"/hr")
print("Background multi rate is "+ str(backgroundMultis/backgroundTime) +"/hr")
