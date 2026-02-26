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

def single_run_v1(path,i):
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
    slowdaq_time_trigger_window = slowdaq_time[compress_idx - 80:compress_idx + time_offset]
    slowsaq_valve_trigger_window = slowdaq_Valve[compress_idx - 80:compress_idx + time_offset]

    return (True, a_fit, t_fit, a_err, t_err)
    

    
    
    
    





































































































































































