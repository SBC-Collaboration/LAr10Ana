from __future__ import division

import copy
import re

import numpy as np
import scipy.signal

def extend_window(w, r):
    # Inputs:
    #   w: An array of 2 elements. Normally, this will be a window like [t1, t2]
    #   r: A float used as a ratio to extend w
    # Outputs: A rescaled version of w
    mp = 0.5*(w[1]+w[0])  # Midpoint
    new_len = (w[1]-w[0])*(1+r)  # Length of new window
    return [mp-new_len/2, mp+new_len/2]

def freq_filter(freqs, lower=None, upper=None):
    # Inputs:
    #   freqs: An array of frequency bins
    #   lower: The lower frequency to cut-off at
    #   upper: The upper frequency to cut-off at
    # Outputs: An array of indices where the the frequency in freqs is between lower and upper
    if lower is None and upper is None:
        return freqs
    if lower is None:
        return np.where([x <= upper for x in freqs])
    if upper is None:
        return np.where([x >= lower for x in freqs])
    return np.where([lower <= x <= upper for x in freqs])


def closest_index(arr, el):
    # Inputs:
    #   arr: A 1-dimensional array
    #   el: Any element
    # Outputs: The FIRST index of the item in arr that is closest to el.
    # Notes: Arr does NOT have to be sorted.
    return np.argmin(np.abs(arr-el))


def spectrum_sums(spectrum, fr, n, lowerf=None, upperf=None):
    # Inputs:
    #   spectrum: The output 2d spectrum from a spectogram
    #   fr: A list of frequency bins corresponding to the spectrum
    #   n: Number of bins
    #   lowerf: The lower frequency to cut-off at
    #   upperf: The upper frequency to cut-off at
    # Outputs: A compressed 1d array where each element is the sum of a bin from spectrum, only counting
    #          frequencies between lowerf and upperf
    out = []
    good_indices = freq_filter(fr, lowerf, upperf)
    for subn in range(n):
        out.append(np.trapezoid(spectrum[good_indices[0], subn], dx=np.mean(np.diff(fr))))
    return out


def rescale_window(w1, w2):
    # Inputs:
    #   w1: An array with 2 elements
    #   w2: An array with 2 elements
    # Outputs: A rescaled version of w2 so tha the endpoints of w2 match w1 but the number of elements remain the same
    y1, y2 = min(w1), max(w1)
    x1, x2 = min(w2), max(w2)
    if x1 == x2:
        return 0*w2
    a = (y1-y2)/(x1-x2)
    b = (x1*y2-x2*y1)/(x1-x2)
    return a*w2+b


def corr_signal(tau, dt, t0, n, fit_type=0, shift=10):
    # Inputs:
    #   tau: Time constant on exponential decay
    #   dt: Step size for the x-axis
    #   t0: Where the exponential signal will start. Not important when used with correlation
    #   N: Number of points requested
    #   fit_type: The type of signal to create. See corr_signal_type_templates.py for a better explanation.
    #               fit_type = 0 --> Exponential decay
    #               fit_type = 1 --> Constant 1 followed by exponential decay (continuous)
    #               fit_type = 2 --> Linear increase followed by exponential decay
    #               fit_type = 3 --> Log increase followed by exponential decay
    #               fit_type = 4 --> 0 value followed by an exponential decrease. Discontinuous.
    # Outputs:
    #   t: t-values for plotting
    #   y: y-values of our filter signal.
    # After careful analysis, we've determined that there reaches a point in the filtered piezo signal that
    # exhibits a sharp increase followed by an exponential decay. This function returns a brief exponential
    # decay function for use with convolution/correlation.
    shift = int(np.ceil(shift))
    t = np.linspace(t0, t0+dt*n, n)
    y = np.exp(-(t-t0)/tau)
    ycopy = copy.deepcopy(y)
    if fit_type == 0:
        pass
    elif fit_type == 1:
        for subn in range(len(y) - shift):
            y[subn+shift] = ycopy[subn]
        y[0:shift] = 1
    elif fit_type == 2:
        for subn in range(len(y) - shift):
            y[subn + shift] = ycopy[subn]
        y[0:shift] = (t[0:shift] - t0)/(shift*dt)
    elif fit_type == 3:
        for subn in range(len(y) - shift):
            y[subn + shift] = ycopy[subn]
            y[0:shift] = np.log((t[0:shift] + 1 - t0)) / np.log(shift*dt + 1)
    elif fit_type == 4:
        for subn in range(len(y) - shift):
            y[subn+shift] = ycopy[subn]
            y[0:shift] = 0
    return t, y

def find_t0_from_corr(corrt, corry):
    # Inputs:
    #   corrt: Time-values of the correlation signal
    #   corry: Y-values of the correlation signal
    # Outputs: The time of the maximum in corry such that corrt is less than or equal to 0.
    # n = np.where(corrt >= 0)
    # corry[n] = 0
    return corrt[np.argmax(corry)]

def spectrogram(piezo_waveform, timebase):
    dt = np.mean(np.diff(timebase))
    return scipy.signal.spectrogram(piezo_waveform, fs=1./dt, nfft=512, noverlap=450,
                                              mode="psd", window="hann", nperseg=512)

def find_peakt0(raw_piezo, times, t0_win, f_low, f_high, n_sample_baseline):
    filtered_piezo = BandPass2(raw_piezo - np.mean(raw_piezo[:n_sample_baseline]), f_low, f_high)
    dt = times[1] - times[0]
    t0_win_ix = np.intp(np.round((t0_win - times[0]) / dt))
    peak_index = np.argmax(np.abs(filtered_piezo[t0_win_ix[0]:t0_win_ix[1]])) + t0_win_ix[0]
    return times[peak_index]

def calculate_t0(piezo_waveform, piezo_timebase, tau, n_sample_baseline=1000,
                 lower=20000, upper=40000, piezo_fit_type=0):
    # Inputs:
    #   piezo_waveform: Piezoelectric waveform
    #   piezo_timebase: The times of each element in the piezo_waveform
    #   tau: The time constant we are trying to fit to the exponential decay that occurs
    #        immediately after the bubble forms
    #   n_sample_baseline: number of samples to use as baseline
    #   lower: The lower frequency threshold for cutting off the spectrogram
    #   upper: The upper frequency threshold for cutting off the spectrogram
    #   piezo_fit_type: The type of fit to use when trying to match the filtered piezo signal. Defaults to 0.
    #   view_plots: Boolean. If true, will display some plots for analysis.
    # Outputs: A dictionary of results for the Acoustic Analysis.
    try:
        timebase = piezo_timebase
        dt = np.mean(np.diff(timebase))
        fr, bn, sp = spectrogram(piezo_waveform, timebase)
        n = len(bn)
        sp_sums = spectrum_sums(sp, fr, n, lower, upper)
        sp_sums = scipy.signal.medfilt(sp_sums)
        textent = [min(timebase), max(timebase)]
        rescaled_t = rescale_window(textent, bn)
        corr_dt = np.mean(np.diff(rescaled_t))
        corr_n = 1000
        corr_t, corr_y = corr_signal(tau, corr_dt, rescaled_t[0], corr_n, fit_type=piezo_fit_type)
        corr = np.correlate(sp_sums, corr_y, "same")
        corr_t = rescaled_t - 0.5 * corr_n * corr_dt
        test_t0 = find_t0_from_corr(corr_t, corr) # This is the t0 we begin to look backwards from

        # Establish a baseline for our lookback algorithm
        #   But first we take the log of the [integrated] spectrogram signal
        log_sp_sums = np.log(sp_sums)

        baseline = np.average(log_sp_sums[:n_sample_baseline])
        baseline_rms = np.std(log_sp_sums[:n_sample_baseline])

        test_t0_index = np.argmin(np.abs(rescaled_t - test_t0))
        print("INITIAL T0", test_t0, test_t0_index)
        rescaled_dt = np.mean(np.diff(rescaled_t))
        t_thresh = 100e-6
        n_lookback = int(np.floor(t_thresh/rescaled_dt))
        pts_lookbacked_sofar = 0

        while True:
            to_test = log_sp_sums[test_t0_index-n_lookback-pts_lookbacked_sofar:test_t0_index-pts_lookbacked_sofar]

            if np.all(to_test<(baseline+5*baseline_rms)):
                break
            pts_lookbacked_sofar += 1
            if test_t0_index-n_lookback-pts_lookbacked_sofar <= 0:
                pts_lookbacked_sofar = -1
                break
        if pts_lookbacked_sofar != -1:
            t0 = rescaled_t[test_t0_index-pts_lookbacked_sofar] + rescaled_dt/2

        else:
            t0 = np.nan
        return t0
    except Exception as e:
        raise
        return np.nan

def BandPass2(yd, f_low, f_high):
    fband = np.array([f_low, f_high])
    b, a = scipy.signal.butter(2, fband / (2.5e6 / 2.0), btype='bandpass', output='ba')
    yd_f = scipy.signal.filtfilt(b, a, yd)
    return yd_f

def CalcPiezoE(yd, td, t_wins, f_bins, t0):
    piezoE = np.zeros((t_wins.shape[0],
                       f_bins.shape[0] - 1),
                      dtype=np.float64) + np.nan

    if np.isnan(t0):
        return piezoE

    dt = td[1] - td[0]
    t_wins = t_wins + t0
    t_wins_ix = np.intp(np.round((t_wins - td[0]) / dt))
    t_wins_ix[t_wins_ix < 0] = 0
    t_wins_ix[t_wins_ix > td.shape[0]] = td.shape[0]


    for i_win in range(t_wins.shape[0]):
        this_yd = yd[t_wins_ix[i_win][0]:t_wins_ix[i_win][1]]
        if len(this_yd) < 2:
            continue
        fft_amp = np.fft.rfft(this_yd)
        fft_pow = (np.abs(fft_amp) ** 2) * dt / len(this_yd)

        df = 1 / (dt * len(this_yd))
        fd = df * (np.arange(len(fft_amp), dtype=np.float64) + 1)
        f_bins_ix = np.intp(np.round((f_bins / df) - 1))
        f_bins_ix[f_bins_ix < 0] = 0
        f_bins_ix[f_bins_ix > len(fft_amp)] = len(fft_amp)

        fft_en = fft_pow * (fd ** 2)
        for i_f in range(len(f_bins) - 1):
            piezoE[i_win, i_f] = df *\
                np.sum(fft_en[f_bins_ix[i_f]:f_bins_ix[i_f + 1]])

    return piezoE



def AcousticAnalysis(ev, tau=0, piezo_fit_type=0,
                     f_high=np.float64(40e3), f_low=np.float64(6e3),
                     t0_win=np.float64([0, 0.1]),
                     n_sample_baseline=np.intp(1e4), 
                     t_wins=np.float64([[-0.1, 0.1]]),
                     f_bins=np.float64([1e2, 1e3, 1e4, 1e5]),
                     corr_lowerf=20000, corr_upperf=40000):
    # Inputs:
    #   ev: Event data (from GetEvent)
    #   tau: The expected time-constant of the exponential decay from the filtered piezo signal
    #   piezo_fit_type: Piezoelectric response fit type
    #   f_high: High frequency used for calculating Piezo Energy
    #   f_low: Low frequency used for calculating Piezo Energy
    #   t0_win: Time window to expect peak piezoelectric response
    #   n_sample_baseline: Number of samples at start of acoustic waveform to use for sampling baseline
    #   t_wins: Time windows for calculating Piezo energy
    #   f_bins: Frequency bins for calculating Piezo energy
    #   corr_lowerf: Lower frequency for cutting off the spectrogram of the Piezo signal
    #   corr_upperf: Higher frequency for cutting off the spectrogram of the Piezo signal
    # Outputs: A dictionary of values for the various variables we are interested in. For a list, take a look
    #          at default_output right below.

    default_output = dict()
    out = default_output

    try:
        if not ev["acoustics"]["loaded"]:
            return default_output

        wvfs = ev["acoustics"]["Waveform"]
        rnge = ev["acoustics"]["Range"]
        dcoffset = ev["acoustics"]["DCOffset"]

        dt = 1/(float(ev["run_control"]["acous"]["sample_rate"][:-5])*1e6) # convert MHz -> Hz
        times = np.arange(wvfs.shape[-1])*dt
  
    except:
        return out

    out["bubble_t0"] = np.zeros(wvfs.shape[1], dtype=np.float64) + np.nan
    out["peak_t0"] = np.zeros(wvfs.shape[1], dtype=np.float64) + np.nan
    out["piezoE"] = np.zeros((wvfs.shape[1], t_wins.shape[0], f_bins.shape[0]-1),
                                 dtype=np.float64) + np.nan

    wvfs = (wvfs.T/rnge.T - dcoffset.T).T # convert to mV, subtract offset

    for i_piezo in range(wvfs.shape[1]):
        try:
            raw_piezo = wvfs[0, i_piezo, :]

            out["peak_t0"][i_piezo] = find_peakt0(raw_piezo, times, t0_win, f_low, f_high, n_sample_baseline)
            t0 = calculate_t0(raw_piezo, times, n_sample_baseline=n_sample_baseline,
                              tau=tau, lower=corr_lowerf, upper=corr_upperf,
                              piezo_fit_type=piezo_fit_type)
            out["bubble_t0"][i_piezo] = t0
            t0_index = closest_index(times, t0)

            out["piezoE"][i_piezo,:,:] = CalcPiezoE(raw_piezo, times, t_wins, f_bins, times[t0_index])

        except Exception as e:
            raise

    return out


if __name__ == "__main__":
    from GetEvent import GetEvent

    TEST_RUN = "/exp/e961/data/SBC-25-daqdata/20250611_1/"
    TEST_EVENT = 0

    tau = 0.0038163479219674467

    out = AcousticAnalysis(GetEvent(TEST_RUN, TEST_EVENT), tau=tau)
    print(out)
