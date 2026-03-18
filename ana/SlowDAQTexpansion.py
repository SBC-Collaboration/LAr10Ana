#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.signal import firwin, filtfilt
from scipy.optimize import least_squares

def piecewise_with_t(params, x):
    # pressure fitting function
    a, t = params
    y = np.where(x < t, 0, a*(x - t))
    return y

def residuals_with_t(params, x, y):
    # chi square
    return piecewise_with_t(params, x) - y


"""
v0 by Ryan
output (Expansion success(bool),  t_start_daq, t_compression , expansion_time)
 expansion success: if the expansion is successful. i.e. pressure acheived at target pressure 1000 ms before the compression
 t_start_daq[ms]: start time of expansion point in slowdaq channel
 t_compression[ms]: compression time. This synchronize the trigger signal in all channels. It is when SERVO valve reaches 100% output
 in slowDAQ channel. It matches 800ms timestamp in acoustic channel with pre-trigger on and 40th camera image timestamp in camera chanle with pre-trigger on
 expansion_time[ms]: the difference between t_start and t_compression. 
 
 
 Notice: 
 1.there is slightly difference between lifetime(timewindow between expansion begins and bubble formation) and expansion time because
 usually bubble formation time is slightly earlier than compression time in order of 100 ms. 
 2. How to get the lifetime: lifetime =  expansion time - (800-t0).
 t0 is the reconstructed bubble formation time in acoustic channel and it is noted as t0_fitting=t_fit in PressureT0.py
 """



def SlowDAQTexpansionFinding(ev, expansion=False,  t_start_daq=0, t_compression = 0, expansion_time = 0):

    default_output = dict(expansion=expansion, t_start_daq=t_start_daq,  t_compression = t_compression, expansion_time = expansion_time
                          )
    try:


        slowdaq_PT = data["slow_daq"]['PT1101']
        slowdaq_Valve = data["slow_daq"]['SERVO3321_OUT']
        slowdaq_time = [i * 10 for i in range(len(slowdaq_PT))]

        # find compression time
        compress_idx = 0
        for i in range(len(data["slow_daq"]['SERVO3321_OUT'])):
            if data["slow_daq"]['SERVO3321_OUT'][i] > 75:
                # print(i, data["slow_daq"]['SERVO3321_OUT'][i], data["slow_daq"]['time_ms'][i])
                t_compression = data["slow_daq"]['time_ms'][i]
                compress_idx = i
                break

        # # find t_start -PT version
        # find t_start
        time_cut = -10
        slowdaq_PT = data["slow_daq"]['PT2121']
        slowdaq_time = data["slow_daq"]['time_ms']
        t_start_time_window = slowdaq_time[:compress_idx + time_cut]
        pressure_cut_compress = slowdaq_PT[:compress_idx + time_cut]
        reverse_time_start = t_start_time_window[::-1]
        reverse_time_start = [int(x) for x in reverse_time_start]
        reverse_pressure = pressure_cut_compress[::-1]
        time_width = reverse_time_start[-2] - reverse_time_start[-1]
        p_set = float(data["event_info"]['pset_hi'])
        if reverse_pressure[0] < p_set + 0.05: # if expansion acheived the target pressure
            for i in range(len(reverse_pressure)):
                if reverse_pressure[i] < p_set + 0.05 and reverse_pressure[i + 1] > p_set + 0.05:
                    interpolation_part = (p_set + 0.05 - reverse_pressure[i]) * (
                                reverse_time_start[i + 1] - reverse_time_start[i]) / (
                                                     reverse_pressure[i + 1] - reverse_pressure[i])
                    interpolation_time = reverse_time_start[i] + interpolation_part
                    # print((p_set+0.05-reverse_pressure[i]),(reverse_time_start[i+1]-reverse_time_start[i]),(reverse_pressure[i+1]-reverse_pressure[i]))
                    # print(interpolation_part, interpolation_time)
                    t_start_daq = interpolation_time

                    break


        else:
            # expansion fails to acheive target pressure
            t_start_daq = t_compression

        expansion_time = t_compression - t_start_daq
        if expansion_time < 1000: # both short expansion or failure to acheive target pressure will trigger the clause
            output = dict(expansion=False, t_start_daq=t_start_daq,  t_compression = t_compression, expansion_time = expansion_time
                          )
        else:
            output = dict(expansion=True, t_start_daq=t_start_daq, t_compression=t_compression,
                          expansion_time=expansion_time
                          )


        return output
    except:
        return default_output


if __name__ =="__main__":
    from sbcbinaryformat import Streamer, Writer
    import numpy as np
    import matplotlib.pyplot as plt

    from GetEvent import GetEvent

    from ana import AcousticT0
    from scipy.signal import firwin, filtfilt
    from scipy.optimize import least_squares
    import importlib

    # data = GetEvent("/exp/e961/app/users/runze/data/20251120_12/", 3,strictMode=False) # success event
    data = GetEvent("/exp/e961/app/users/runze/data/20251120_12/", 5, strictMode=False)  # success event
    result  = SlowDAQTexpansionFinding(data)
    print(result)