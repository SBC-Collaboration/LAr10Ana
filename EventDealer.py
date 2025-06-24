import os
import re
import shutil
import time

import numpy as np
import copy
import numpy.matlib

from ana.EventAnalysis import EventAnalysis as eva
from AcousticT0 import AcousticAnalysis as aa
from ana.ExposureAnalysis import ExposureAnalysis as expa 
from ana.SiPMPulses import SiPMPulses as sa

from GetEvent import GetEvent
from sbcbinaryformat import Streamer, Writer

ANALYSES = {
    "event": eva,
    # "acoustic": aa,
    "exposure": expa,
    "scintillation": sa,
}

def BuildEventList(rundir, first_event=0, last_event=-1):
    # Inputs:
    #   rundir: Directory for the run
    #   first_event: Index of first event
    #   last_event: Index of last_event
    # Outputs: A sorted list of events from rundir
    # Possible Issues: If we ever move to a different naming scheme, this needs to be re-written.
    eventdirlist = []
    for f in os.listdir(rundir):
        if f.isdigit() and os.path.isdir(os.path.join(rundir, f)):
            eventdirlist.append(int(f))
    eventdirlist = np.array(sorted(eventdirlist))
    if last_event >= 0:
        eventdirlist = eventdirlist[eventdirlist <= last_event]
    return eventdirlist

# map numpy dtype.str names to those expected by sbcbinaryformat
def dname(s):
    s = s[1:] # remove leading carat
    if s == "f4": s = "f" # set float name
    if s == "f8": s = "d" # set double name

    return s

def ProcessSingleRun(rundir, dataset='SBC-25', recondir='.', process_list=None):
    # Inputs:
    #   rundir: Location of raw data
    #   dataset: Indicator used for filtering which analyses to run
    #   recondir: Location of recon data/where we want to output our binary files
    #   process_list: List of analyses modules to run. example: ["acoustic", "event", ""]
    # Outputs: Nothing. Saves binary files to recondir.
    if process_list is None:
        process_list = []  # This is needed since lists are mutable objects. If you have a default argument
                           # as a mutable object, then the default argument can *change* across multiple
                           # function calls since the argument is created ONCE when the function is defined.
    runname = os.path.basename(rundir)
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
    eventlist = BuildEventList(rundir)

    for ev in eventlist:
        t0 = time.time()
        print('Starting event ' + runname + '/' + str(ev))

        data = GetEvent(rundir, ev)
        print('Time to load event:  '.rjust(35) +
              str(time.time() - t0) + ' seconds')
        npev = np.array([ev], dtype=np.int32)

        for p in process_list:
            t1 = time.time()
            out[p].append(ANALYSES[p](data, **parameter_config[p]))
            out[p][-1]['runid'] = runid
            out[p][-1]['ev'] = npev
            et = time.time() - t1
            print(('%s analysis:  ' % p).rjust(35) + str(et) + ' seconds')

        print('*** Full event analysis ***  '.rjust(35) +
              str(time.time() - t0) + ' seconds')

    # save everything
    for p in process_list:
        column_names = list(out[p][0].keys())
        dtypes = [dname(out[p][0][c].dtype.str) for c in column_names]
        sizes = [list(out[p][0][c].shape) for c in column_names]
        writer = Writer(os.path.join(run_recondir, p + runname + ".bin"), column_names, dtypes, sizes)
        for evind in range(len(out[p])):
            writer.write(dict([(c, out[p][evind][c]) for c in column_names]))
    return
