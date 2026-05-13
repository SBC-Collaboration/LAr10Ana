from sbcbinaryformat import Streamer, Writer
import numpy as np
import matplotlib.pyplot as plt

from GetEvent import GetEvent
from GetEvent import GetScint

from scipy import signal

import csv
import pandas as pd

import os

from sklearn.linear_model import RANSACRegressor
from sklearn.linear_model import LinearRegression



def _unwrap_caen_timestamp(ts, max_ts=2**31):
    ts = np.asarray(ts, dtype=np.int64)

    # Detect rollovers
    rollovers = np.diff(ts, axis=-1, prepend=0) < 0

    # Cumulative count of rollovers
    rollover_count = np.cumsum(rollovers, axis=-1)
    return ts + rollover_count * max_ts


def find_offset_correlation(t1_original, t2_original, bin_width=0.00001):
    t1 = np.asarray(t1_original).copy()
    t2 = np.asarray(t2_original).copy()
    t1_0 = t1[0]
    t2_0 = t2[0]
    t1 -= t1_0
    t2 -= t2_0
    
    t_min = min(t1[0], t2[0])
    t_max = max(t1[-1], t2[-1])
    bins = np.arange(t_min, t_max, bin_width)
    hist1, _ = np.histogram(t1, bins=bins)
    hist2, _ = np.histogram(t2, bins=bins)
    
    correlation = signal.correlate(hist1, hist2, mode='full')
    lag = signal.correlation_lags(len(hist1), len(hist2), mode='full')
    lag = lag*bin_width + t1_0 - t2_0
    
    best_lag = lag[np.argmax(correlation)]
    
    return best_lag, correlation, lag




def scint_t0(data):
    results = {}

    results['Failed'] = False
    results['latch_time_corrected'] = np.nan 
    results['pT0_in_scint_time'] = np.nan 
    results['biggest_pulse_pt0_20ms'] = np.nan 
    results['idx_biggest_pulse_pt0_20ms'] = np.nan 
    results['scint_time_biggest_pulse_pt0_20ms'] = np.nan 
    results['biggest_pulse_random_20ms'] = np.nan 
    results['idx_biggest_pulse_random_20ms'] = np.nan 
    results['scint_time_biggest_pulse_random_20ms'] = np.nan 
    results['biggest_pulse_pt0_40ms'] = np.nan 
    results['idx_biggest_pulse_pt0_40ms'] = np.nan 
    results['scint_time_biggest_pulse_pt0_40ms'] = np.nan 
    results['biggest_pulse_random_40ms'] = np.nan 
    results['idx_biggest_pulse_random_40ms'] = np.nan 
    results['scint_time_biggest_pulse_random_40ms'] = np.nan 

    
    ### Digiscope
    
    digi_data = data['digiscope']
    
    di_array = digi_data['DI']
    first_num = digi_data['DI'][:,0]

    bit0 = ((first_num >> 0) & 1).astype(bool)

    rising_edges = np.where((~bit0[:-1]) & (bit0[1:]))[0] + 1

    if len(rising_edges) == 0:
        results['Failed'] = 'No digiscope latch found'
        return results
        
    latch_idx = rising_edges[-1]
    bit1 = ((first_num >> 1) & 1).astype(bool)
    bit2 = ((first_num >> 2) & 1).astype(bool)

    caen_trig_idx = np.where(bit2)[0]
    if len(caen_trig_idx) == 0:
        results['Failed'] = 'No digiscope triggers found'
        return results


    unwrapped_timestamps_digi1 = _unwrap_caen_timestamp(digi_data['t_ticks'])
    unwrapped_timestamps_digi = _unwrap_caen_timestamp(unwrapped_timestamps_digi1)

    latch_timestamp = unwrapped_timestamps_digi[latch_idx]
    latch_time = latch_timestamp / 40e6

    first_caen_timestamp = unwrapped_timestamps_digi[caen_trig_idx[0]]
    last_caen_timestamp = unwrapped_timestamps_digi[caen_trig_idx[-1]]

    latch_minus_first_caen = latch_timestamp-first_caen_timestamp
    latch_time_after_first_caen = latch_minus_first_caen / 40e6
    last_caen_minus_latch_s = (last_caen_timestamp-latch_timestamp) / 40e6

    caen_livetime = (last_caen_timestamp - first_caen_timestamp) / 40e6

    timestamps_digi = unwrapped_timestamps_digi[caen_trig_idx]
    digi_tsec = (timestamps_digi / 40e6)





    ### Scintillation
    
    scint_data = data['scintillation']

    sample_rate = 125e6
    
    TriggerTimeTag = _unwrap_caen_timestamp(scint_data['TriggerTimeTag'])

    if len(TriggerTimeTag) == 0:
        results['Failed'] = 'No scint triggers found'
        return results
        
    scint_livetime = (TriggerTimeTag[-1]-TriggerTimeTag[0]) / sample_rate
    scint_tsec = ((TriggerTimeTag-TriggerTimeTag[0]) / sample_rate)





    ### Digi/Scint Correlation
    
    t, corr, lag = find_offset_correlation(digi_tsec, scint_tsec, bin_width=0.001) # Initial rough correlation, before clock drift fix

    digi_aligned = np.array(digi_tsec) - t
    scint = np.array(scint_tsec)

    
    # Nearest neighbor, residuals
    
    digi_aligned.sort()
    scint.sort()

    idx = np.searchsorted(scint, digi_aligned)

    idx_left = np.clip(idx - 1, 0, len(scint) - 1)
    idx_right = np.clip(idx, 0, len(scint) - 1)

    diff_left = digi_aligned - scint[idx_left]
    diff_right = digi_aligned - scint[idx_right]

    choose_left = np.abs(diff_left) < np.abs(diff_right)
    residuals = np.where(choose_left, diff_left, diff_right)

    worst_idx = np.argmax(np.abs(residuals))
    worst_residual = residuals[worst_idx]

    threshold = 0.002  # seconds, initial outlier cutoff
    residuals_filtered = residuals[np.abs(residuals) < threshold]


    mask = np.abs(residuals) <= threshold

    digi_aligned_filtered = digi_aligned[mask]
    residuals_filtered = residuals[mask]


    # Fit clock drift
    X = digi_aligned_filtered.reshape(-1, 1)
    y = residuals_filtered
    if len(X) < 2:
        results['Failed'] = 'Not enough points for RANSAC fit'
        return results
        
    ransac = RANSACRegressor(LinearRegression(),
                             residual_threshold=0.0000005)  # 0.5us

    ransac.fit(X, y)

    inlier_mask = ransac.inlier_mask_
    outlier_mask = ~inlier_mask

    m = ransac.estimator_.coef_[0]  # clock drift slope
    b = ransac.estimator_.intercept_


    x_in = digi_aligned_filtered[inlier_mask]
    y_in = residuals_filtered[inlier_mask]

    resid_to_fit = y_in - (m * x_in + b)

    resid_us = resid_to_fit * 1e6



    # Recheck nearest neighbor, residuals

    digi_inliers = digi_aligned_filtered[inlier_mask]
    residuals_inliers = residuals_filtered[inlier_mask]

    digi_corrected = (1 - m) * digi_inliers - b

    idx = np.searchsorted(scint, digi_corrected)

    idx_left = np.clip(idx - 1, 0, len(scint) - 1)
    idx_right = np.clip(idx, 0, len(scint) - 1)

    diff_left = digi_corrected - scint[idx_left]
    diff_right = digi_corrected - scint[idx_right]

    choose_left = np.abs(diff_left) < np.abs(diff_right)
    residuals_final = np.where(choose_left, diff_left, diff_right)

    

    latch_time_corrected = ((1 - m) * (latch_time-t) - b) * 1000 # ms

    results['latch_time_corrected'] = latch_time_corrected

    
    
    ### PressureT0

    if "pressure_t0" in data["analysis"]:
        pressureT0_results = data["analysis"]["pressure_t0"]

        pressureT0 = pressureT0_results['t0_fitting']  # ms
        
        if pressureT0 == 0:
            results['Failed'] = 'pT0 = 0'
            return results
            
        pt0_rel_to_trig = -1*((data['run_control']['acous']['pre_trig_len'] * 1000) - pressureT0)
    
    else:
        results['Failed'] = 'pT0 results not found'
        return results

    pressureT0_in_corrected_time = latch_time_corrected + pt0_rel_to_trig # ms
    results['pT0_in_scint_time'] = pressureT0_in_corrected_time



    ### Scint hits
    
    if "scintillation" in data["analysis"]:
        scintillation_hits = data["analysis"]["scintillation"]

    else:
        results['Failed'] = 'Scintillation hits results not found'
        return results


    hit_area = np.asarray(scintillation_hits['hit_area'], dtype=float)
    hit_sums = np.nansum(hit_area, axis=1)





    # Biggest pulse within 20ms of pT0

    pt0_window = 0.020  # seconds
    
    window_mask = np.abs(scint_tsec - (pressureT0_in_corrected_time/1000)) <= pt0_window
    
    window_indices = np.where(window_mask)[0]
    
    if len(window_indices) == 0:
        biggest_pulse_within_20ms_pt0 = np.nan
        index_biggest_pulse_within_20ms_pt0 = np.nan
        scint_time_biggest_pulse_within_20ms_pt0 = np.nan
    else:
        window_values = hit_sums[window_indices]
    
        local_max_idx = np.nanargmax(window_values)
        index_biggest_pulse_within_20ms_pt0 = window_indices[local_max_idx]
    
        biggest_pulse_within_20ms_pt0 = window_values[local_max_idx]

        scint_time_biggest_pulse_within_20ms_pt0 = scint_tsec[index_biggest_pulse_within_20ms_pt0]
    

    results['biggest_pulse_pt0_20ms'] = biggest_pulse_within_20ms_pt0
    results['idx_biggest_pulse_pt0_20ms'] = index_biggest_pulse_within_20ms_pt0
    results['scint_time_biggest_pulse_pt0_20ms'] = scint_time_biggest_pulse_within_20ms_pt0 * 1000





    # Biggest pulse within random +-20ms window (excluding pT0 window)

    pt0_center = pressureT0_in_corrected_time / 1000
    
    # Random window only before pT0 - 20ms
    random_center_min = pt0_window
    random_center_max = pt0_center - (2* pt0_window)
    
    biggest_pulse_random_20ms = np.nan
    index_biggest_pulse_random_20ms = np.nan
    scint_time_biggest_pulse_random_20ms = np.nan
    
    if random_center_max > random_center_min: # failed pT0 gives pT0 = 0
    
        random_center = np.random.uniform(random_center_min,random_center_max)
        random_window_mask = np.abs(scint_tsec - random_center) <= pt0_window
        random_window_indices = np.where(random_window_mask)[0]

    
        if len(random_window_indices) > 0:
    
            random_window_values = hit_sums[random_window_indices]
    
            local_max_idx = np.nanargmax(random_window_values)
            index_biggest_pulse_random_20ms = random_window_indices[local_max_idx]
    
            biggest_pulse_random_20ms = random_window_values[local_max_idx]
    
            scint_time_biggest_pulse_random_20ms = scint_tsec[index_biggest_pulse_random_20ms]
    
    
    results['biggest_pulse_random_20ms'] = biggest_pulse_random_20ms
    results['idx_biggest_pulse_random_20ms'] = index_biggest_pulse_random_20ms
    results['scint_time_biggest_pulse_random_20ms'] = scint_time_biggest_pulse_random_20ms * 1000





    # Biggest pulse within 40ms of pT0

    pt0_window = 0.040
    
    window_mask = np.abs(scint_tsec - (pressureT0_in_corrected_time/1000)) <= pt0_window
    
    window_indices = np.where(window_mask)[0]
    
    if len(window_indices) == 0:
        biggest_pulse_within_40ms_pt0 = np.nan
        index_biggest_pulse_within_40ms_pt0 = np.nan
        scint_time_biggest_pulse_within_40ms_pt0 = np.nan
    else:
        window_values = hit_sums[window_indices]
    
        local_max_idx = np.nanargmax(window_values)
        index_biggest_pulse_within_40ms_pt0 = window_indices[local_max_idx]
    
        biggest_pulse_within_40ms_pt0 = window_values[local_max_idx]

        scint_time_biggest_pulse_within_40ms_pt0 = scint_tsec[index_biggest_pulse_within_40ms_pt0]
    

    results['biggest_pulse_pt0_40ms'] = biggest_pulse_within_40ms_pt0
    results['idx_biggest_pulse_pt0_40ms'] = index_biggest_pulse_within_40ms_pt0
    results['scint_time_biggest_pulse_pt0_40ms'] = scint_time_biggest_pulse_within_40ms_pt0 * 1000





    # Biggest pulse within random +-40ms window (excluding pT0 window)

    pt0_center = pressureT0_in_corrected_time / 1000
    
    # Random window only before pT0 - 40ms
    random_center_min = pt0_window
    random_center_max = pt0_center - (2* pt0_window)
    
    biggest_pulse_random_40ms = np.nan
    index_biggest_pulse_random_40ms = np.nan
    scint_time_biggest_pulse_random_40ms = np.nan
    
    if random_center_max > random_center_min: # failed pT0 gives pT0 = 0
    
        random_center = np.random.uniform(random_center_min,random_center_max)
        random_window_mask = np.abs(scint_tsec - random_center) <= pt0_window
        random_window_indices = np.where(random_window_mask)[0]
    
    
        if len(random_window_indices) > 0:
    
            random_window_values = hit_sums[random_window_indices]
    
            local_max_idx = np.nanargmax(random_window_values)
            index_biggest_pulse_random_40ms = random_window_indices[local_max_idx]
    
            biggest_pulse_random_40ms = random_window_values[local_max_idx]
    
            scint_time_biggest_pulse_random_40ms = scint_tsec[index_biggest_pulse_random_40ms]
    
    
    results['biggest_pulse_random_40ms'] = biggest_pulse_random_40ms
    results['idx_biggest_pulse_random_40ms'] = index_biggest_pulse_random_40ms
    results['scint_time_biggest_pulse_random_40ms'] = scint_time_biggest_pulse_random_40ms * 1000
    
    
    
    return results
    

    
    

    
    