# This should be run after successfully running convert_raw_to_npy_run_by_run.py
# This merges the npy files for each raw run into a single npy file called raw_events.npy

# Usage example: python merge_raw_run_npy.py /bluearc/storage/30l-16-data
# Produces npy files for each run for navigation from raw data, to be used by PED event display
# may need to source /coupp/data/home/coupp/PEDsvn/setup_ped_paths.sh

from glob import glob
import numpy as np
import os
import re
import time
import sys

skip = ['timestamp', 'livetime', 'piezo_max(3)', 'piezo_min(3)', 'piezo_starttime(3)', 'piezo_endtime(3)', 'piezo_freq_binedges(9)', 'acoustic_neutron', 'acoustic_alpha', 'scanner_array(2)', 'scan_source_array(2)', 'scan_nbub_array(2)', 'scan_trigger_array(2)', 'scan_comment_array(2)', 'scaler(8)', 'led_max_amp(8)', 'led_max_time(8)', 'null_max_amp(8)', 'first_hit(8)', 'last_hit(8)', 'max_amps(8)', 'max_times(8)', 'nearest_amps(8)', 'nearest_times(8)', 'numtrigs(8)', 'numpretrigs(8)', 'scan_comment_array(2)']
dtypes = {'s': 'U12', 'd': 'i4', 'f': 'f4', 'e': 'f4'}  # map fscanf format to numpy datatype

if len(sys.argv) != 2:
    print('Should be 1 argument.')
    print('Usage: python convert.py <dir-where-raw-data-is>')
    exit()

# print('Number of arguments:' + str(len(sys.argv)) + 'arguments.')
# print('Argument List:' + str(sys.argv))

raw_directory = str(sys.argv[1])
# print('raw_directory = ' + raw_directory)

def natural_sort(things):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(things, key=alphanum_key)


print('Starting now')
start = time.time()
runs = [os.path.basename(x) for x in glob(os.path.join(raw_directory, "20*.npy"))]
runs = natural_sort(runs)

all_events = np.array([], dtype=[('run', 'U12'), ('ev', 'i4'), ('reco index', 'i4')])
counter = 0
for run in runs:
    print(run)
    try:
        run_events = np.load(os.path.join(raw_directory, run))
        all_events = np.concatenate([all_events, run_events])
        print("Added " + str(len(run_events)) + " events from " + run)
        counter = counter + 1
    except Exception as e:
        print(e)
        print("Failed to add events from " + run)

print("Saving " + str(len(all_events)) + " events from " + str(counter) + " runs to raw_events.npy")
np.save(os.path.join(raw_directory, 'raw_events'), all_events)        

print('finished in {:.0f} seconds'.format(time.time() - start))



