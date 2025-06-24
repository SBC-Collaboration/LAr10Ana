# This creates reco_events.npy from merged_all.txt and re-indexes the raw_events.npy file for fast cuts on reco variables in the EventDisplay
# Note this also creates reco_events_all.npy which contains all the reco events, while reco_events.npy is culled if necessary to only include events present in the raw data
# (this should mostly only be relevant to those with local copies of a smaller fraction of the raw datasets)

from glob import glob
import numpy as np
import os
import re
import time
import sys

skip = ['timestamp', 'livetime', 'piezo_max(3)', 'piezo_min(3)', 'piezo_starttime(3)', 'piezo_endtime(3)', 'piezo_freq_binedges(9)', 'acoustic_neutron', 'acoustic_alpha', 'scanner_array(2)', 'scan_source_array(2)', 'scan_nbub_array(2)', 'scan_trigger_array(2)', 'scan_comment_array(2)', 'scaler(8)', 'led_max_amp(8)', 'led_max_time(8)', 'null_max_amp(8)', 'first_hit(8)', 'last_hit(8)', 'max_amps(8)', 'max_times(8)', 'nearest_amps(8)', 'nearest_times(8)', 'numtrigs(8)', 'numpretrigs(8)', 'scan_comment_array(2)']
dtypes = {'s': 'U12', 'd': 'i4', 'f': 'f4', 'e': 'f4'}  # map fscanf format to numpy datatype

if len(sys.argv) < 3:
    print('Should be 2 to 4 arguments.')
    print('Usage: python convert.py <dir-where-reco-data-is> <dir-where-npy-are> <optional-arg-merged-filename> <optional-arg-user-date>')
    exit()

reco_directory = str(sys.argv[1])
npy_location = str(sys.argv[2])
print("npy files will be put at: " + npy_location)
merged_filename = 'merged_all.txt'
if len(sys.argv) > 3:
    merged_filename = str(sys.argv[3])    

print('Using merged filename: ', merged_filename)

user_date = ''
if len(sys.argv) == 5:
    user_date = str(sys.argv[4])
    
def natural_sort(things):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(things, key=alphanum_key)

def load_reco(filename):
    path = os.path.join(reco_directory, filename)

    with open(path) as f:
        files = f.readline()
        fields = f.readline().split()
        formats = f.readline().split()
        
        dt = []
        columns = []
        column = 0
        for field in fields:
            format = formats[column]
            for x in dtypes:
                if x in format:
                    dtype = dtypes[x]
            match = re.search('.*\((.*)\)', field)
            if match:
                dimensions = [int(x) for x in match.groups()[0].split(',')]
                length = np.prod(dimensions)
                types = formats[column:column + length]
                if len(set(types)) > 1:
                    message = 'cannot parse {} with types {} because mixed types are not supported'
                    raise NotImplementedError(message.format(field, types))
                if (len(dimensions) == 1) and (field not in skip):  # skip loading multidimensional arrays to save memory
                    columns.extend(list(range(column, column + length)))
                    dt.append((field, (dtype, dimensions)))
                column += length
            else:
                if field not in skip:
                    dt.append((field, dtype))
                    columns.append(column)
                column += 1

    skip_header = 6
    new_events = np.genfromtxt(path, dtype=np.dtype(dt), skip_header=skip_header, usecols=columns)
    
    return new_events

def load_raw(filename, reco_all):
    try:
        raw = np.load(os.path.join(npy_location, filename))
        print("Saving backup of raw_events.npy as raw_events_bkp.npy")
        np.save(os.path.join(npy_location, 'raw_events_bkp'), raw)
    except FileNotFoundError:
        print('Cannot find raw data!')
        #raw = np.array([], dtype=[('run', 'U12'), ('ev', 'i4'), ('reco index', 'i4')])
        return None

    raw_events = []
    reco_events = []
    reco_index = 0
    n_reco_evt = 0
    # print(raw)
    # print("reco_all.type: ", type(reco_all))
    # print("reco_all,shape: ", reco_all.shape)
    # print("reco_all,dtype: ", reco_all.dtype)
    for row in raw:
        run = row['run']
        event = row['ev'] 
        # print(run, event)
        matches = np.argwhere((reco_all['run'] == run) & (reco_all['ev'] == int(event))).flatten()
        # print( 'len matches: ', len(matches))
        # print( 'reco_index: ', reco_index)
        if len(matches) > 0:
            raw_events.append((run, event, reco_index))
            reco_events.extend( reco_all[matches] )
            # print('reco_all[matches]: ', reco_all[matches])
            # print('added: ', reco_events[-1])
            # reco_index = reco_index + 1
            reco_index = reco_index + len(matches)
            n_reco_evt += 1
        else:
            raw_events.append((run, event, -1))

        
    # print("reco_events.type: ", type(reco_events))
    # print(reco_events)
    new_raw = np.array(raw_events, dtype=[('run', 'U12'), ('ev', 'i4'), ('reco index', 'i4')])
    new_reco = np.array(reco_events)#, dtype=reco_all[0].dtype)
    # print(new_reco)
    # print("new_reco.type: ", type(new_reco))
    # print("new_reco,shape: ", new_reco.shape)
    # print("new_reco,dtype: ", new_reco.dtype)
    # diff_events = np.setdiff1d(new_raw, raw)
    print('Number of events found in read raw npy file: {}'.format(len(raw)))
    print('Number of events to be put in reindexed raw npy file: {}'.format(len(new_raw)))
    print('Number of lines found in read reco npy file: {}'.format(len(reco_all)))
    print('Number of lines to be put in matched reco npy file: {}'.format(len(new_reco)))
    print('Number of events to be put in matched reco npy file: {}'.format(n_reco_evt))
    # print('Length of diff between old and new: {}'.format(len(diff_events)))
    
    return new_raw, new_reco


print('starting now')
start = time.time()

if len(sys.argv) < 5:
    try:
        reco_all = load_reco(merged_filename)
        print("Saving the full reco events npy file as reco_events_all.npy")
        np.save(os.path.join(npy_location, 'reco_events_all'), reco_all)
        try:
            raw, reco = load_raw('raw_events.npy', reco_all)
            print("Saving the re-indexed raw_events.npy")
            np.save(os.path.join(npy_location, 'raw_events'), raw)
            print("Saving the reco that is culled-to-match the raw as reco_events.npy")
            np.save(os.path.join(npy_location, 'reco_events'), reco)
        except Exception as e:
            print('Cannot find raw or reco data, skipping creation of raw_events.npy and reco_events.npy')
            print(e)
    except Exception as e:
        print('Problem with merged_all.txt, Aborting.')
        print(e)
else:
    try:
        reco_all = load_reco(merged_filename)
        print("Saving the full reco events npy file as reco_events_all_{}.npy".format(user_date))
        np.save(os.path.join(npy_location, 'reco_events_all_{}'.format(user_date)), reco_all)
        try:
            raw, reco = load_raw('raw_events.npy', reco_all)
            print("Saving the re-indexed raw_events.npy")
            np.save(os.path.join(npy_location, 'raw_events'), raw)
            print("Saving the reco that is culled-to-match the raw as reco_events_{}.npy".format(user_date))
            np.save(os.path.join(npy_location, 'reco_events_{}'.format(user_date)), reco)
        except Exception as e:
            print('Cannot find raw or reco data, skipping creation of raw_events.npy and reco_events.npy')
            print(e)
    except Exception as e:
        print('Problem with merged_all.txt, Aborting.')
        print(e)

print('finished in {:.0f} seconds'.format(time.time() - start))
