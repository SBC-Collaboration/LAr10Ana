#!/usr/bin/env python3

import numpy as np

PTs = [2121, 1101]

def exposure(PT, time, bins):
    centers = (bins[1:] + bins[:-1])/2
    binned,_ = np.histogram(PT, bins=bins)
    imode = np.argmax(binned)
    DELTA = 0
    pressure = np.mean(PT[PT < centers[imode+DELTA]])
    livetime = np.sum(np.diff(time[PT < centers[imode+DELTA]]))

    return pressure, livetime/1e3

def ExposureAnalysis(ev, pressure_bins=np.linspace(0, 10, 201), PTs=PTs):
    output = {}
    for PT in PTs:
        output["PT%i_pressure" % PT] = 0.
        output["PT%i_livetime" % PT] = 0.

    #try:
    if ev is None or not (ev['event_info']['loaded'] and ev['slow_daq']['loaded']):
        return output

    dt = np.mean(np.diff(ev['slow_daq']['time_ms']))

    for PT in PTs:
        P, T = exposure(ev["slow_daq"]["PT%i" % PT], ev['slow_daq']['time_ms'], pressure_bins)
        output["PT%i_pressure" % PT] = P
        output["PT%i_livetime" % PT] = T

    return output
