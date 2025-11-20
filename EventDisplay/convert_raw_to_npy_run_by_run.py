# This version creates the raw npy files for each run separately, for each that does not already exist
# This is a 2020 redesign for running on Compute Canada resources in which most runs will be in tar format

# Usage examples:
#    python convert.py /bluearc/storage/30l-16-data
#    python convert.py /bluearc/storage/30l-16-data /bluearc/storage/30l-16-data/npy
# Produces npy files for navigation from raw data, to be used by PED event display
# may need to source /coupp/data/home/coupp/PEDsvn/setup_ped_paths.sh

from glob import glob
import numpy as np
import os
import re
import time
import sys
import tarfile
import shutil
import getpass
import zipfile

tar_postfix = '.tar'
# tar_postfix = '.tar.gz'
# tar_postfix = '.tgz'
tar_postfix_len = len(tar_postfix)

delete_untar = True

if len(sys.argv) < 2 or len(sys.argv) > 3:
    print('Should be 1 or 2 arguments.')
    print('To put npy files in same dir as raw data: python convert.py <dir-where-raw-data-is>')
    print('To specify npy file dir: python convert.py <dir-where-raw-data-is> <dir-where-npy-will-be-put>')
    exit()

# print('Number of arguments:' + str(len(sys.argv)) + 'arguments.')
# print('Argument List:' + str(sys.argv))

raw_directory = str(sys.argv[1])
# extraction_path = raw_directory
extraction_path = '../EventDisplay/scratch'
extraction_path = os.path.join(extraction_path, getpass.getuser())
#Will make scratch diretory if doesn't exists
if not os.path.exists(extraction_path):
    os.makedirs(extraction_path)
print("Tarfile extraction path will be: " + extraction_path)

npy_location = raw_directory
if len(sys.argv) == 3:
    npy_location = str(sys.argv[2])
print("npy files will be put at: " + npy_location)

warning_list = []
new_npy_file_list = []
new_untar_list = []
new_zipread_list = []

def natural_sort(things):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(things, key=alphanum_key)

def validate(events, run_folder_path):
    res = []
    for run, event, index in events:
        path = os.path.join(run_folder_path, str(event))
        if os.path.isdir(path):
            res.append((run, event, index))
        else:
            print('  WARNING: dir not found at {}'.format(path))

    return np.array(res, dtype=events.dtype)

def make_npy_of_run(run, run_folder_path):
    events = []
    for event in natural_sort(glob(os.path.join(run_folder_path, '[0-9]*/'))):
        event = os.path.basename(event.strip(os.sep))
        events.append((run, event, -1)) # note we are setting the reco index to -1 here because there is no reco at this point

    print('  Events in run {}: {}'.format(run,len(events)))
    events = np.array(events, dtype=[('run', 'U12'), ('ev', 'i4'), ('reco index', 'i4')])
    events = validate(events, run_folder_path)

    if events.size == 0:
        print(f"  No events found in run {run}; skipping npy generation.")
        return

    try:
        np.save(os.path.join(npy_location, run), events)
        new_npy_file_list.append(run)
    except:
        warning = "WARNING: failed to produce npy file for " + run
        print("  " + warning)
        warning_list.append(warning)

def validate_zip(events, run_folder_path):
    archive = zipfile.ZipFile(run_folder_path, 'r')

    res = []
    for run, event, index in events:
        path = os.path.join(run, str(event))
        try:
            file = archive.open(path, 'r')
            res.append((run, event, index))
        except:
            print('  WARNING: dir not found at {}'.format(path))
        
    return np.array(res, dtype=events.dtype)
        
def make_npy_of_zip_run(run, run_folder_path):
    archive = zipfile.ZipFile(run_folder_path, 'r')

    basenames = []
    for dir in archive.namelist():
        if dir.endswith('/'):
            basename = str(os.path.basename(os.path.normpath(str(dir))))
            if basename.isnumeric():
                basenames.append(basename)
    sortednames = natural_sort(basenames)
    events = []
    for event in sortednames:
        events.append((run, event, -1)) # note we are setting the reco index to -1 here because there is no reco at this point

    print('  Events in run {}: {}'.format(run,len(events)))
    events = np.array(events, dtype=[('run', 'U12'), ('ev', 'i4'), ('reco index', 'i4')])
    events = validate_zip(events, run_folder_path)

    if events.size == 0:
        print(f"  No events found in zipped run {run}; skipping npy generation.")
        return

    try:
        np.save(os.path.join(npy_location, run), events)
        new_npy_file_list.append(run)
    except:
        warning = "WARNING: failed to produce npy file for " + run
        print("  " + warning)
        warning_list.append(warning)
        
        
print('Starting now')
start = time.time()
all_runs = [os.path.basename(x) for x in glob(os.path.join(raw_directory, "202*"))]
all_runs_npy = [os.path.basename(x) for x in glob(os.path.join(npy_location, "202*"))]

# Filter to only include years >= 2025 for SBC
runs = [x for x in all_runs if x[:4].isdigit() and int(x[:4]) >= 2025]
runs_npy = [x for x in all_runs_npy if x[:4].isdigit() and int(x[:4]) >= 2025]

for run in runs:
    path = os.path.join(raw_directory, run)
    print("Checking " + run)
    if run.endswith(".npy"):
        print("It's an npy file")
        if run[:-4] in runs:
            print(" Found run " + run[:-4] + " corresponding to npy file")
        elif run[:-4] + tar_postfix in runs:
            print(" Found run " + run[:-4] + tar_postfix + " corresponding to npy file")
        else:
            warning = "WARNING: npy file found with no corresponding run directory or tar"
            print(" " + warning)
            warning_list.append(warning)
    elif run.endswith(tar_postfix):
        print("It's a tar file")
        run_nopostfix = run[:-tar_postfix_len]
        if run_nopostfix in runs:
            print(" Found run folder " + run_nopostfix + " corresponding to tar file. No need to untar")
        elif run_nopostfix + ".npy" in runs:
            print(" Found npy file " + run_nopostfix + ".npy corresponding to tar file. No need to untar")
        elif run_nopostfix + ".npy" in runs_npy:
            print(" The npy file for this run already exists in npy dir. Ignoring.")
        else:
            print(" Did not find run folder " + run_nopostfix + " or npy file corresponding to tar file. Untarring...")
            try:
                t = tarfile.open(path, 'r')
                t.extractall(extraction_path)
                run = run_nopostfix
                print(" Untar successful. Proceeding to generate npy file.")
                new_untar_list.append(run)
                run_folder_path = os.path.join(extraction_path, run)
                # Make npy file if no npy is found and had to untar
                make_npy_of_run(run, run_folder_path)
                os.system('chmod -R 777 ' + run_folder_path)
                shutil.rmtree(run_folder_path)
            except Exception as e:
                warning = "WARNING: Untar of file " + path + " failed"
                print(" " + warning)
                print(e)
                warning_list.append(warning)
    elif run.endswith('.zip'):
        print("It's a zip file")
        run_nopostfix = run[:-4]
        if run_nopostfix in runs:
            print(" Found run folder " + run_nopostfix + " corresponding to tar file. No need to untar")
        elif run_nopostfix + ".npy" in runs:
            print(" Found npy file " + run_nopostfix + ".npy corresponding to tar file. No need to untar")
        elif run_nopostfix + ".npy" in runs_npy:
            print(" The npy file for this run already exists in npy dir. Ignoring.")
        else:
            print(" Did not find run folder " + run_nopostfix + " or npy file corresponding to zip file. Attempting to get inside...")
            run = run_nopostfix
            try:
                make_npy_of_zip_run(run, path)
                new_zipread_list.append(run)
            except Exception as e:
                warning = "WARNING: zip file access " + path + " failed"
                print(" " + warning)
                print(e)
                warning_list.append(warning)                
    elif os.path.isdir(path):
        print("It's a run directory")
        if run + ".npy" in runs:
            print(" The npy file for this run already exists in raw dir. Ignoring.")
        elif run + ".npy" in runs_npy:
            print(" The npy file for this run already exists in npy dir. Ignoring.")
        else:
            print(" No npy file found: generating npy file.")
            run_folder_path = os.path.join(raw_directory, run)
            # Make npy file if no npy is found and didn't have to untar 
            make_npy_of_run(run, run_folder_path)
    else:
        warning = "WARNING: I don't know what this file is. Ignoring " + run
        print(warning)
        warning_list.append(warning)

        
print("#################################")
print('finished in {:.0f} seconds'.format(time.time() - start))
print("####### List of warnings: #######")
for warning in warning_list:
    print(warning)
print("####### List of untars: #######")
for untarred_run in new_untar_list:
    print(untarred_run)
print("####### List of read zips: #######")
for zip_run in new_zipread_list:
    print(zip_run)
print("####### List of npy files made: #######")
for new_npy_file in new_npy_file_list:
    print(new_npy_file)