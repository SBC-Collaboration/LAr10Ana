import os
import re
import shutil
import time
import sys
import numpy as np
import copy
import numpy.matlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ana.EventAnalysis import EventAnalysis as eva
from ana.AcousticT0 import AcousticAnalysis as aa
from ana.ExposureAnalysis import ExposureAnalysis as expa 
from ana.SiPMPulses import SiPMPulsesBatched as sa
from ana.ScintRate import ScintillationRateAnalysis as sra

from GetEvent import GetEvent, NEvent
from sbcbinaryformat import Streamer, Writer

ANALYSES = {
    "event": eva,
    # "acoustic": aa,
    "exposure": expa,
    "scintillation": sa,
    "scint_rate": sra,
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

    # output data and defaults
    out = dict([(p, []) for p in process_list])
    defaults = dict([(p, ANALYSES[p](None)) for p in process_list])

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

    for ev in eventlist:
        t0 = time.time()
        print('Starting event ' + runname + '/' + str(ev))

        data = GetEvent(rundir, ev, strictMode=False)
        print('Time to load event:  '.rjust(35) + f"{time.time()-t0:.6f} seconds")
        npev = np.array([ev], dtype=np.int32)

        for p in process_list:
            t1 = time.time()

            try:
                out[p].append(ANALYSES[p](data, **parameter_config[p]))
            except Exception as e:
                print("Analysis %s failed on event %i with error: %s" % (p, ev, str(e)))
                out[p].append(defaults[p])

            out[p][-1]['runid'] = runid
            out[p][-1]['ev'] = npev
            et = time.time() - t1
            print(('%s analysis:  ' % p).rjust(35) + f"{et:.6f} seconds")

        print('*** Full event analysis ***  '.rjust(35) + f"{time.time()-t0:.6f} seconds\n")

    # save everything
    for p in process_list:
        if not out[p]:
            continue
        
        column_names = list(out[p][0].keys())
        dtypes = []
        for c in column_names:
            val = out[p][0][c]
            # Convert to numpy array if it isn't already
            if not isinstance(val, np.ndarray):
                val = np.array(val)
            dtypes.append(dname(val.dtype.str))

        # squeeze sizes
        sizes = [list(np.squeeze(out[p][0][c]).shape) for c in column_names]
        # set default
        sizes = [s if len(s) else [1] for s in sizes]
        # for outputs with a sub-event number, fix the sizes
        sizes = [s[1:] if len(s) > 1 else s for s in sizes]

        writer = Writer(os.path.join(run_recondir, f"{p}.sbc"), column_names, dtypes, sizes)
        for evind in range(len(out[p])):
            writer.write(dict([(c, np.squeeze(out[p][evind][c])) for c in column_names]))
    return
