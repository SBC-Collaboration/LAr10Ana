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
v1 by Ryan
output (t0_fiting, a_fitting, t0_uncertainty, a_uncertainty, chi_square value)
 t0_fiting[ms], a_fitting[bara/ms**2], t0_uncertainty[ms], a_uncertainty[bara/ms**2], chi_square value: in the fitting, the pressure change is fit to flat constant 0 + linear increasing function
 described by the t0: the junction point and a, the linear slope. The t0 here is in the timescale of acoustic channel
 """
def PressureT0Finding(ev, t0_fitting = 0, a_fitting=0, t0_sigma = 0, a_sigma = 0, t0_chi_sq=0):

    default_output = dict(t0_fitting = t0_fitting, a_fitting=a_fitting,
                          t0_sigma = t0_sigma, a_sigma = a_sigma, t0_chi_sq=t0_chi_sq)
    try:


        wvfs = ev["acoustics"]["Waveforms"]

        wvfs_psi = wvfs * (-35 / 2 ** 16)

        piezo0 = wvfs_psi[0, 7, :]
        xlimit = [0, 800]
        ylimit = [-22.5, -20]
        ylimit = [-20, -18]
        # first check if the expansion success
        average_window_expansion = 100


        average_window = 50  # every 50 data, do the average
        n_chunked = (len(piezo0) // average_window) * average_window  # chunk data
        piezo0 = piezo0[:n_chunked]

        # time in miliseconds
        total_time = len(piezo0)
        time_list_ms = [i / 1e3 for i in range(0, total_time, 1)]
        time_list_ms = np.array(time_list_ms[:n_chunked])

        # max time to be fit
        time_max_raw = max(time_list_ms)


        piezo0 = piezo0.reshape(-1, average_window).mean(axis=1)
        time_list_ms = time_list_ms.reshape(-1, average_window).mean(axis=1)

        # print("time length", len(time_list_ms))
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
                fitting_ending_indx = int(min(i - int(10000 / average_window), int(time_max_raw*1000 / average_window)))
                break
            else:
                fitting_ending_indx = int(time_max_raw*1000 / average_window)
        # print("index", ending_indx, time_list_ms[ending_indx])
        pressure_before_fit = piezo0_filtered[starting_indx:ending_indx]
        time_range = time_list_ms[starting_indx:ending_indx]
        # print("fitting indx", starting_indx, fitting_ending_indx,ending_indx)

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
                # print(index)
                break
        # 8e5

        # fit with function
        # f =c when x<t0,
        # f=a(x-t0)**2+c when x>t0

        slope_before_fit = slope0_filtered[starting_indx:fitting_ending_indx]
        time_fitting_range = time_list_ms[starting_indx:fitting_ending_indx]
        # print("node0")
        # initial guesses:
        a0 = 8e-6
        t0 = time_list_ms[fitting_ending_indx]-100  # t0 initial guess around 600 ms

        # c0 = np.mean(pressure_before_fit[:100]) # first 100 data average
        p0 = [a0, t0]

        # optionally set bounds: a>=0 , t within x-range, c in pressure range
        bounds = ([0, min(time_fitting_range)], [np.inf, max(time_fitting_range)])
        # print(bounds)

        res = least_squares(residuals_with_t, p0, args=(time_fitting_range, slope_before_fit), bounds=bounds)
        a_fit, t_fit = res.x
        # print("fitted a, t:", a_fit, t_fit)

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

        output = dict(t0_fitting=t_fit, a_fitting=a_fit,
                      t0_sigma=t_err, a_sigma=a_err, t0_chi_sq=sigma2)
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

    # data = GetEvent("/exp/e961/data/SBC-25-unpacked/20251113_9/", 3,strictMode=False) # success event
    data = GetEvent("/exp/e961/data/SBC-25-unpacked/20251112_18/", 33, strictMode=False)  # failure event
    # data = GetEvent("/exp/e961/data/SBC-25-unpacked/20251114_0/", 1, strictMode=False)

    result  = PressureT0Finding(data)
    print(result)
    bad_result = ['20251108_1', '20251114_0', '20251116_2', '20251117_14', '20251123_2', '20260106_3']

    # data_list = ["20251107_41", "20251107_42", "20251108_1", "20251108_4", "20251108_5", "20251110_0", "20251110_1", "20251110_2", "20251111_2", "20251111_3", "20251111_4", "20251111_5", "20251111_6", "20251111_7", "20251111_8", "20251112_18", "20251112_19", "20251113_0", "20251113_1", "20251113_9", "20251113_10", "20251113_11", "20251114_0", "20251114_1", "20251114_6", "20251114_36", "20251114_37", "20251115_0", "20251115_1", "20251115_2", "20251115_3", "20251115_4", "20251115_5", "20251116_1", "20251116_2", "20251117_0", "20251117_1", "20251117_12", "20251117_14", "20251117_15", "20251118_0", "20251118_1", "20251118_2", "20251118_3", "20251118_4", "20251118_8", "20251118_9", "20251119_0", "20251119_1", "20251119_2", "20251119_3", "20251119_10", "20251120_0", "20251120_12", "20251120_13", "20251121_1", "20251121_5", "20251122_0", "20251122_1", "20251122_2", "20251123_0", "20251123_2", "20251123_3", "20251123_(6", "20251124_(0", "20251125_6", "20251125_7", "20251125_8", "20251126_2", "20251126_3", "20251126_4", "20251126_5", "20251126_7", "20251126_8", "20251127_0", "20251128_0", "20251129_0", "20251130_0", "20251201_3", "20251205_0", "20251206_0", "20251207_0", "20251208_0", "20251210_2", "20251211_0", "20251211_5", "20251211_6", "20251211_7", "20251211_9", "20251211_10", "20251211_11", "20251212_0", "20251212_1", "20251213_6", "20251213_7", "20251215_2", "20251215_3", "20251216_0", "20251216_1", "20251216_13", "20251216_14", "20251216_15", "20251217_0", "20251218_0", "20251218_6", "20251219_0", "20260106_2", "20260106_3", "20260107_3", "20260108_0", "20260108_4", "20260109_0", "20260109_2", "20260109_5", "20260109_10", "20260110_0", "20260110_3", "20260112_2", "20260112_3", "20260113_0", "20260113_4", "20260115_5", "20260116_0", "20260117_0", "20260118_0", "20260119_0", "20260120_0", "20260120_6", "20260120_9", "20260121_1", "20260121_2", "20260122_0", "20260122_3", "20260123_0", "20260123_8", "20260124_0", "20260124_3", "20260125_0", "20260130_2", "20260131_0", "20260131_8", "20260201_0", "20260202_0", "20260202_13", "20260202_16", "20260203_0", "20260203_8", "20260203_10", "20260204_0", "20260205_0", "20260205_5", "20260205_12", "20260206_0", "20260206_8", "20260206_19", "20260206_21", "20260207_0", "20260208_0", "20260209_0", "20260209_5", "20260209_11", "20260210_0", "20260212_0", "20260212_1", "20260212_4", "20260213_1", "20260213_4", "20260214_0", "20260215_0", "20260216_0", "20260217_0", "20260217_7", "20260217_10", "20260218_0", "20260218_9", "20260218_15", "20260219_0","20260220_1", "20260220_5","20260221_0", "20260222_0","20260223_0","20260224_0", "20260225_0","20260225_6", "20260225_11","20260226_0", "20260227_0","20260227_6", "20260227_10","20260228_0","20260228_4", "20260228_7","20260301_0","20260301_3","20260301_6", "20260302_0","20260303_0", "20260303_2","20260303_4"]# test module:
    # bad_list = []
    # for path in data_list:
    #     try:
    #         print(path)
    #         read_path = "/exp/e961/data/SBC-25-unpacked/" + path + "/"
    #         data = GetEvent(read_path, 1, strictMode=False)
    #         result = PressureT0Finding(data)
    #         result_list = list(result.values())
    #         if not any(result_list):
    #             print("bad path", path)
    #             bad_list.append(path)
    #         else:
    #             print(result_list)
    #     except Exception as e:
    #         print("error",e)
    #         continue
    #
    # print(bad_list)
