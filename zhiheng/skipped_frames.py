import numpy as np
import pandas as pd
from pathlib import Path

data_path = Path('/exp/e961/data/SBC-25-daqdata/')
run = '20250827_1'
run_path = data_path/run
events = [p for p in run_path.iterdir() if p.is_dir()]
events = sorted(events, key=lambda p: int(p.name))
cams = ['cam1', 'cam2', 'cam3']
total = {'cam1': 0, 'cam2': 0, 'cam3': 0}
# skipped counts the number of frames taken that have more than 1 skipped frame associated
# gaps counts the total number of skipped frames. So if a frame taken has 2 skipped frames before (30ms gap), then it would add 2 to gaps, and 1 to skipped
skipped = {'cam1': 0, 'cam2': 0, 'cam3': 0}
gaps = {'cam1': 0, 'cam2': 0, 'cam3': 0}
print(f"Run: {run}, {len(events)} events")

for ev in events:
    for cam in cams:
        try:
            info = pd.read_csv(ev/f'{cam}-info.csv')
            total[cam] += len(info)
            frames = info['skipped'][info['skipped'] > 0]
            gaps[cam] += frames.sum()            
            skipped[cam] += frames.count()
            if frames.count() > 0:
                # print event, camera, number of skipped frames in the event, and index, amount skipped per frame.
                print(f"Event {ev.name}, {cam}: {frames.count()} skipped frames: "
                      f"{list(zip(frames.index, map(int, frames.values)))}")
        except Exception as e:
            continue

print("")
for c in range(3):
    cam = cams[c]
    print(f"Camera {c+1} - Total: {total[cam]},\t Skipped: {skipped[cam]} "
          f"({skipped[cam] / total[cam] * 100 if total[cam] > 0 else 0 :.2f}%),\t"
          f"Gaps: {gaps[cam]}")
