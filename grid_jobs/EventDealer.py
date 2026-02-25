# EventDealer.py

import os
import re
import shutil
import time
import sys
import numpy as np
import copy
import numpy.matlib
import gc

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ana.EventAnalysis import EventAnalysis as eva
from ana.AcousticT0 import AcousticAnalysis as aa
from ana.ExposureAnalysis import ExposureAnalysis as expa 
from ana.SiPMPulses import SiPMPulsesBatched as sa
from ana.ScintRate import ScintillationRateBatched as sra
from ana.BubbleFinder import BubbleFinder as bf

from GetEvent import GetEvent, NEvent
from sbcbinaryformat import Streamer, Writer

ANALYSES = {
    "event": eva,
    # "acoustic": aa,
    "exposure": expa,
    "scintillation": sa,
    "scint_rate": sra,
    "bubble": bf
}

def BuildEventList(rundir, maxevt=-1):
    # Inputs:
    #   rundir: Directory for the run
    #   maxevt: Maximum number of events to process
    # Outputs: A sorted list of events from rundir
    # Possible Issues: If we ever move to a different naming scheme, this needs to be re-written.
    eventdirlist = np.array(range(NEvent(rundir)))
    if maxevt >= 0:
        eventdirlist = eventdirlist[eventdirlist < maxevt]
    return eventdirlist

# map numpy dtype.str names to those expected by sbcbinaryformat
def dname(s):
    s = s[1:] # remove leading carat
    if s == "f4": s = "f" # set float name
    if s == "f8": s = "d" # set double name

    return s

def ProcessSingleRun(rundir, dataset='SBC-25', recondir='.', process_list=None, maxevt=-1):
    # Inputs:
    #   rundir: Location of raw data
    #   dataset: Indicator used for filtering which analyses to run
    #   recondir: Location of recon data/where we want to output our binary files
    #   process_list: List of analyses modules to run. example: ["acoustic", "event", ""]
    #   maxevt: Maximum number of events to process
    # Outputs: Nothing. Saves binary files to recondir.
    if process_list is None:
        process_list = []  # This is needed since lists are mutable objects. If you have a default argument
                           # as a mutable object, then the default argument can *change* across multiple
                           # function calls since the argument is created ONCE when the function is defined.
    runname = os.path.basename(rundir).split(".")[0]
    runid_str = runname.split('_')
    runid = np.int32(runid_str)
    # run_recondir = os.path.join(recondir, runname)
    run_recondir = recondir

    # Make output directory
    if not os.path.isdir(run_recondir):
        os.mkdir(run_recondir)

    # get processes
    process_list = [p.lower().strip() for p in process_list]
    
    # process parameter configuraton
    parameter_config = dict([(p, {}) for p in process_list])

    # Acoustic analysis parameters
    if dataset == "SBC-25":
        if "acoustic" in process_list:
            # G.P. -- COPIED FROM SBC-17. TODO: revisit at some point
            tau_peak = 0.0025884277467056165  # <-- This comes from TauResultAnalysis.py (J.G.)
            tau_average = 0.0038163479219674467  # <-- This also ^^
            # lower_f = 20000 OLD...
            # upper_f = 40000
            lower_f = 1000
            upper_f = 25000
            piezo_fit_type = 0
            
            parameter_config["acoustic"]["tau"] = tau_average
            parameter_config["acoustic"]["piezo_fit_type"] = piezo_fit_type
            parameter_config["acoustic"]["corr_lowerf"] = lower_f
            parameter_config["acoustic"]["corr_upperf"] = upper_f
            
    print("Starting run " + rundir)
    eventlist = BuildEventList(rundir, maxevt=maxevt)

    # Create writers before event loop
    writers = {}

    for ev in eventlist:
        t0 = time.time()
        print('Starting event ' + runname + '/' + str(ev))

        try:
            data = GetEvent(rundir, ev, strictMode=False)
        except Exception as e:
            print(f"Failed to load event {ev} with error: {e}. Skipping event.")
            continue

        print('Time to load event:  '.rjust(35) + f"{time.time()-t0:.6f} seconds")
        npev = np.array([ev], dtype=np.int32)

        for p in process_list:
            t1 = time.time()

            # Skip analysis if data not loaded
            if (p == "scint_rate" or p == "scintillation") and not data["scintillation"]["loaded"]:
                print(f"Skipping {p} analysis -- scintillation data not loaded.")
                continue
            elif p == "exposure" and not (data["event_info"]["loaded"] and data["slow_daq"]["loaded"]):
                print(f"Skipping {p} analysis -- event info data not loaded.")
                continue
            elif p == "acoustic" and not data["acoustic"]["loaded"]:
                print(f"Skipping {p} analysis -- acoustic data not loaded.")
                continue
            elif p == "event" and not data["event_info"]["loaded"]:
                print(f"Skipping {p} analysis -- event info data not loaded.")
                continue
            elif p == "bubble" and not data["cam"]["loaded"]:
                print(f"Skipping {p} analysis -- event info data not loaded.")
                continue

            try:
                result = ANALYSES[p](data, **parameter_config[p])
            except Exception as e:
                print("Analysis %s failed on event %i with error: %s" % (p, ev, str(e)))
                continue
            result['runid'] = runid
            result['ev'] = npev
            
            # create writer if it doesn't exist
            if p not in writers:
                column_names = list(result.keys())
                dtypes = []
                sizes = []
                
                for c in column_names:
                    val = result[c]
                    if not isinstance(val, np.ndarray):
                        val = np.array(val)
                    dtypes.append(dname(val.dtype.str))
                    
                    if p == "scint_rate" or p == "bubble":
                        shape = list(np.atleast_1d(val).shape)
                    else:
                        shape = list(np.squeeze(val).shape)
                    shape = shape if len(shape) else [1]
                    shape = shape[1:] if len(shape) > 1 else shape
                    sizes.append(shape)
                
                writers[p] = Writer(os.path.join(run_recondir, f"{p}.sbc"), column_names, dtypes, sizes)

            # Write to file
            column_names = list(result.keys())
            writers[p].write(dict([(c, np.squeeze(result[c])) for c in column_names]))
            del result
            
            et = time.time() - t1
            print(('%s analysis:  ' % p).rjust(35) + f"{et:.6f} seconds")
        
        del data
        gc.collect()

        print('*** Full event analysis ***  '.rjust(35) + f"{time.time()-t0:.6f} seconds\n")

    # delete all writers
    for p in process_list:
        if p in writers:
            del writers[p]
    
    return

if __name__ == "__main__":
    if len(sys.argv) > 1:
        ProcessSingleRun(
            rundir=sys.argv[1],
            recondir=sys.argv[2],
            process_list = ["event", "exposure", "scintillation", "scint_rate", "bubble"])
    else:
        ProcessSingleRun(
            rundir="/exp/e961/data/SBC-25-daqdata/20260221_0.tar",
            recondir="/home/zsheng/test", # Use your own directory for testing~
            process_list = ["event", "bubble"])