from sbcbinaryformat import Streamer, Writer
import numpy as np
import matplotlib.pyplot as plt

from GetEvent import GetEvent

from ana import AcousticT0 
from scipy.signal import firwin, filtfilt
from scipy.optimize import least_squares
import importlib

importlib.reload(AcousticT0)


def piecewise_with_t(params, x):
    a, t = params
    y = np.where(x < t, 0, a*(x - t))
    return y

def residuals_with_t(params, x, y):
    return piecewise_with_t(params, x) - y

def single_run_daq_v1(path,i):
    # outpuut(bool if expansion success, fit of slope value, fit of t0 value, uncerntainty of a value, uncerntainty of t0 value)
    TEST_RUN = "/exp/e961/data/SBC-25-daqdata/20250611_1/"
    TEST_RUN2 = "/exp/e961/app/users/runze/data/20251120_12/"
    # TEST_RUN2 = "/exp/e961/data/users/gputnam/SBC-25-daqdata-test/20251103_1/"
    TEST_EVT = i
    
    data = GetEvent(path, i,strictMode=False)

    slowdaq_PT = data["slow_daq"]['PT1101']
    slowdaq_Valve = data["slow_daq"]['SERVO3321_OUT']
    slowdaq_time = [i * 10 for i in range(len(slowdaq_PT))]

    # find compression time
    compress_idx = 0
    for i in range(len(data["slow_daq"]['SERVO3321_OUT'])):
        if data["slow_daq"]['SERVO3321_OUT'][i] > 75:
            # print(i, data["slow_daq"]['SERVO3321_OUT'][i], data["slow_daq"]['time_ms'][i])
            compress_time = data["slow_daq"]['time_ms'][i]
            compress_idx = i
            break

    

     # # find t_start -PT version
    # find t_start
    time_cut  = -10
    slowdaq_PT = data["slow_daq"]['PT2121']
    slowdaq_time = data["slow_daq"]['time_ms']
    t_start_time_window = slowdaq_time[:compress_idx+time_cut]
    pressure_cut_compress = slowdaq_PT[:compress_idx+time_cut]
    reverse_time_start = t_start_time_window[::-1]
    reverse_time_start = [int(x) for x in reverse_time_start]
    reverse_pressure =pressure_cut_compress[::-1]
    time_width = reverse_time_start[-2]-reverse_time_start[-1]
    p_set = float(data["event_info"]['pset_hi'])
    # print("set", p_set, "mode",p_mode)
    if reverse_pressure[0]<p_set+0.05:
        for i in range(len(reverse_pressure)):
            if reverse_pressure[i]<p_set+0.05 and reverse_pressure[i+1]>p_set+0.05:
                interpolation_part = (p_set+0.05-reverse_pressure[i])*(reverse_time_start[i+1]-reverse_time_start[i])/(reverse_pressure[i+1]-reverse_pressure[i])
                interpolation_time = reverse_time_start[i] + interpolation_part
                # print((p_set+0.05-reverse_pressure[i]),(reverse_time_start[i+1]-reverse_time_start[i]),(reverse_pressure[i+1]-reverse_pressure[i]))
                # print(interpolation_part, interpolation_time)
                t_start = interpolation_time
        
        
                break
       
                
    else:
        t_start = compress_time

    time_diff = compress_time-t_start
    if time_diff< 850:
        return (False,  t_start, compress_time, compress_time ,0)
    
    # find t0 in valve channel find the turning point before compress time

    startpoint = 10000
    endpoint = 10750
    time_offset = 0
    pressure_diff = slowdaq_PT[:-1] - slowdaq_PT[1:]

    # filter servo valve
    numtaps = 10  # filter length (longer = sharper cutoff)
    Fs = 1000  # sampling rate
    Fc = 10  # low pass filter in Hz
    fir = firwin(numtaps, Fc, window='hamming', fs=Fs)

    valve_filtered = filtfilt(fir, [1.0], slowdaq_Valve[compress_idx - 80:compress_idx + time_offset])
    # backwards 80 points and 10ms per points so 800 ms
    slowdaq_time_trigger_window = slowdaq_time[compress_idx - 80:compress_idx + time_offset]
    slowsaq_valve_trigger_window = slowdaq_Valve[compress_idx - 80:compress_idx + time_offset]

    time_window= slowdaq_time[compress_idx-80:compress_idx+time_offset]
    time_window_shifted = []
    for i in range(len(time_window)):
        time_window_shifted.append(time_window[i]- time_window[0])

    valve_delta = valve_filtered[1:]-valve_filtered[:-1]
    reverse_time = time_window_shifted[::-1]
    reverse_time = [int(x) for x in reverse_time]
    reverse_valve_delta =valve_delta[::-1]
    time_width = reverse_time[-2]-reverse_time[-1]
    for i in range(len(reverse_valve_delta)):
        if reverse_valve_delta[i]<0 and reverse_valve_delta[i+1]>0:
           
            interpolation_part = (0-reverse_valve_delta[i])*(reverse_time[i+1]-reverse_time[i])/(reverse_valve_delta[i+1]-reverse_valve_delta[i])
            interpolation_time = reverse_time[i] + interpolation_part
            t0_valve_acoustic = interpolation_time
            
            break
    t0_valve_daq = t0_valve_acoustic+ time_window[0]

   
    return (True,  t_start, compress_time, t0_valve_daq ,t0_valve_acoustic)

    
    

    
    
    
    





































































































































































