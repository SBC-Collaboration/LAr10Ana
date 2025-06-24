#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 14:59:54 2024

@author: cdahl
"""

import numpy as np

PTs = [3308, 3309, 3311, 3314, 3320, 3332, 3333]

def ExposureAnalysis(ev,
                     PressureEdge=13.5+np.linspace(0,87,88, dtype=np.float64),
                     TimeEdge=np.float64([0, 5, 10, 20, 30]),
                     BubbleTime=np.float64([-0.205, -0.15])
                     ):

    EventPressure=-1+np.zeros((len(PTs),len(BubbleTime)), dtype=np.float64)
    Exposure=-1+np.zeros((len(PTs),len(PressureEdge)-1,len(TimeEdge)+1), dtype=np.float64)

    default_output = dict(PressureEdge=PressureEdge,
                          TimeEdge=TimeEdge,
                          BubbleTime=BubbleTime,
                          Exposure=Exposure,
                          EventPressure=EventPressure
                          )
    #try:
    if ev is None or not (ev['event_info']['loaded'] and ev['slow_daq']['loaded']):
        return default_output

    dt = np.mean(np.diff(ev['slow_daq']['time_ms']))
    # trig_ix = np.nonzero(np.diff(ev['slow_daq']['TriggerLatch'])==1)[0][0]+1
    # exp_ixarray = np.nonzero(np.diff(ev['slow_daq']['TriggerLatch'])==-1)[0] #[-1]+1
    trig_ix = 0 # G.P. TODO: how to look this up
    exp_ix = 0 # ditto

    bubix = np.intp(trig_ix + np.round(BubbleTime/dt))
    lt_end = ev['slow_daq']['time_ms'][trig_ix]-0.5*dt
    lt_start = lt_end - ev['event_info']['ev_livetime']
    ev_start = ev['slow_daq']['time_ms'][0]
    ev_end = ev['slow_daq']['time_ms'][-1]

    time_histedges = np.concatenate((ev_start[None],
                                     # TimeEdge+lt_start,
                                     ev_end[None]))

    for i_p, PT in enumerate(PTs):
        P = ev['slow_daq']["PT%i" % PT]
        EventPressure[i_p] = P[bubix]
        P[:exp_ix] = -99
        P[trig_ix:] = -99
        
        Exposure[i_p] = dt * np.histogram2d(P, ev['slow_daq']['time_ms'],
                                            [PressureEdge, time_histedges])[0]

    output = dict(PressureEdge=PressureEdge,
                  TimeEdge=TimeEdge,
                  BubbleTime=BubbleTime,
                  Exposure=Exposure,
                  EventPressure=EventPressure
                  )
    return output
