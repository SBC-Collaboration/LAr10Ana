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
    # chi sqaure
    return piecewise_with_t(params, x) - y


"""output (Expansion success(bool), t0_fiting, a_fitting, t0_uncertainty, a_uncertainty, chi_square value, t_start_in_daq,t0_in_daq)
 expansion success: if the expansion is successful. i.e. if the pressure achieves stable state before the following bubble formation trigger
 t0_fiting, a_fitting, t0_uncertainty, a_uncertainty, chi_square value: in the fitting, the pressure change is fit to flat constant 0 + linear increasing function
 described by the t0: the junction point and a, the linear slope. The t0 here is in the timescale of acoustic channel
 t_start_in_daq: start time of expansion point in slowdaq channel
 t0_in_daq, t0 value mapping to slowdaq channel"""
def PressureT0Finding(ev, expansion=False, t0_fitting = 0, a_fitting=0, t0_sigma = 0, a_sigma = 0, t0_chi_sq=0, t_start_daq=0, t0_daq = 0):

    default_output = dict(expansion=expansion, t0_fitting = t0_fitting, a_fitting=a_fitting,
                          t0_sigma = t0_sigma, a_sigma = a_sigma, t0_chi_sq=t0_chi_sq, t_start_daq=t_start_daq, t0_daq = t0_daq
                          )
    # try:
    if ev is None or not (ev["acoustics"]["Waveforms"]):
        return default_output

    wvfs = ev["acoustics"]["Waveforms"]

    wvfs_psi = wvfs * (-35 / 2 ** 16)

    piezo0 = wvfs_psi[0, 7, :]
    xlimit = [0, 800]
    ylimit = [-22.5, -20]
    ylimit = [-20, -18]
    # first check if the expansion success
    average_window_expansion = 100
    start_pressure = np.mean(piezo0[:average_window_expansion])
    end_pressure = np.mean(piezo0[600000:600000 + average_window_expansion])
    if start_pressure - end_pressure > 0.3:
        return default_output

    average_window = 50  # every 50 data, do the average
    n_chunked = (len(piezo0) // average_window) * average_window  # chunk data
    piezo0 = piezo0[:n_chunked]

    # time in miliseconds
    total_time = len(piezo0)
    time_list_ms = [i / 1e3 for i in range(0, total_time, 1)]
    time_list_ms = np.array(time_list_ms[:n_chunked])

    piezo0 = piezo0.reshape(-1, average_window).mean(axis=1)
    time_list_ms = time_list_ms.reshape(-1, average_window).mean(axis=1)

    print("time length", len(time_list_ms))
    # add low pass filter
    # assuming 1 microsecond time resolution 1e6Hz
    numtaps = 1000  # filter length (longer = sharper cutoff)
    Fs = 1000000 / average_window  # sampling rate
    Fc = 10  # low pass filter in Hz
    fir = firwin(numtaps, Fc, window='hamming', fs=Fs)

    piezo0_filtered = filtfilt(fir, [1.0], piezo0)
    piezoslope0 = (+piezo0_filtered[1:] - piezo0_filtered[:-1]) * 1e3  # in bar/ms


    starting_indx = int(10000 / average_window)
    # starting 10ms after data collection
    # print(piezoslope0[starting_indx])
    p_rate_range = abs(max(piezoslope0[starting_indx:int(40 * starting_indx)]) - min(
        piezoslope0[starting_indx:int(40 * starting_indx)]))

    ending_indx = 0
    fitting_ending_indx = 0
    hardcut_threshold = p_rate_range * 2
    for i in range(starting_indx, len(piezoslope0), 1):
        if piezoslope0[i] > 2 * hardcut_threshold:
            ending_indx = i
            fitting_ending_indx = i - int(100000)  # modify this
            fitting_ending_indx = int(min(i - int(10000 / average_window), 800000 / average_window))
            break
    # print("index", ending_indx, time_list_ms[ending_indx])
    pressure_before_fit = piezo0_filtered[starting_indx:ending_indx]
    time_range = time_list_ms[starting_indx:ending_indx]
    # print("fitting indx", starting_indx, fitting_ending_indx)

    # add low pass filter
    # assuming 1 microsecond time resolution 1e6Hz
    numtaps = 1000  # filter length (longer = sharper cutoff)
    Fs = 1000000 / average_window  # sampling rate
    Fc = 10  # low pass filter in Hz
    fir = firwin(numtaps, Fc, window='hamming', fs=Fs)

    slope0_filtered = filtfilt(fir, [1.0], piezoslope0)

    # manually choose starting point and end point
    for index in range(len(time_list_ms)):
        if time_list_ms[index] > 800:
            print(index)
            break
    # 8e5

    # fit with function
    # f =c when x<t0,
    # f=a(x-t0)**2+c when x>t0
    slope_before_fit = slope0_filtered[starting_indx:fitting_ending_indx]
    time_fitting_range = time_list_ms[starting_indx:fitting_ending_indx]
    print(time_fitting_range)
    # initial guesses:
    a0 = 8e-6
    t0 = 600  # t0 initial guess around 100ms
    # c0 = np.mean(pressure_before_fit[:100]) # first 100 data average
    p0 = [a0, t0]

    # optionally set bounds: a>=0 , t within x-range, c in pressure range
    bounds = ([0, min(time_fitting_range)], [np.inf, max(time_fitting_range)])
    print(bounds)

    res = least_squares(residuals_with_t, p0, args=(time_fitting_range, slope_before_fit), bounds=bounds)
    a_fit, t_fit = res.x
    print("fitted a, t:", a_fit, t_fit)

    n = len(res.fun)
    p = len(res.x)

    # residual variance (reduced chi^2 estimate)
    sigma2 = np.sum(res.fun ** 2) / (n - p)

    # covariance matrix
    J = res.jac
    cov = sigma2 * np.linalg.inv(J.T @ J)

    # 1-sigma uncertainties
    param_errors = np.sqrt(np.diag(cov))

    a_err, t_err = param_errors

    output = dict(expansion=True, t0_fitting=t_fit, a_fitting=a_fit,
                  t0_sigma=t_err, a_sigma=a_err, t0_chi_sq=sigma2, t_start_daq=0, t0_daq=0
                  )
    return output


if __name__ =="__main__":
    from sbcbinaryformat import Streamer, Writer
    import numpy as np
    import matplotlib.pyplot as plt

    from GetEvent import GetEvent

    from ana import AcousticT0
    from scipy.signal import firwin, filtfilt
    from scipy.optimize import least_squares
    import importlib

    data = GetEvent("/exp/e961/app/users/runze/data/20251120_12/", 3,strictMode=False)
    result  = PressureT0Finding(data)
    print(result)