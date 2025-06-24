# This version creates the raw npy files for a run given as input, and is meant to be run by the analysis chain

# Usage example: python convert.py /bluearc/storage/30l-16-data/run
# Produces npy files for navigation from raw data, to be used by PED event display

from glob import glob
import numpy as np
import os
import re
import time
import sys


if len(sys.argv) != 2:
    print('Should be 1 argument.')
    print('Usage: python convert.py <run-folder-path>')
    exit()

# print('Number of arguments:' + str(len(sys.argv)) + 'arguments.')
# print('Argument List:' + str(sys.argv))

runpath = str(sys.argv[1])
               
def natural_sort(things):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(things, key=alphanum_key)

def validate(events):
    res = []
    for run, event, index in events:
        path = os.path.join(runpath, str(event), 'Event.txt')
        if os.path.isfile(path):
            res.append((run, event, index))
        else:
            print('  WARNING: Event.txt not found at {}'.format(path))

    return np.array(res, dtype=events.dtype)

def make_npy_of_run():
    run = os.path.basename(runpath)
    events = []
    for event in natural_sort(glob(os.path.join(runpath, '[0-9]*/'))):
        event = os.path.basename(event.strip(os.sep))
        events.append((run, event, -1)) # note we are setting the reco index to -1 here because there is no reco at this point

    print('  Events in run {}: {}'.format(run,len(events)))
    events = np.array(events, dtype=[('run', 'U12'), ('ev', 'i4'), ('reco index', 'i4')])
    events = validate(events)
    try:
        np.save(runpath, events)
    except:
        print("WARNING: failed to produce npy file for " + runpath)

        
print('Starting convert_raw_to_npy_for_run.py')
start = time.time()
if os.path.isdir(runpath):
    print("Opening directory " + runpath)
    make_npy_of_run()
else:
    print("WARNING: " + runpath + " does not seem to be a directory. Aborting")
        
print('finished in {:.0f} seconds'.format(time.time() - start))
