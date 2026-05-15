from sbcbinaryformat import Streamer, Writer
import numpy as np
import matplotlib.pyplot as plt

from GetEvent import GetEvent
from GetEvent import GetScint



#def single_expansion_check_run(path,i):
def single_expansion_check_run(data):
    
    wvfs = data["acoustics"]["Waveforms"]
    wvfs_psi = wvfs*(-35/2**16)

    piezo0 = wvfs_psi[0, 7, :]

    average_window = 100
    start_pressure = np.mean(piezo0[:average_window])
    end_pressure = np.mean(piezo0[600000:600000+average_window])
    if start_pressure-end_pressure<0.5:
        return True
    else:
        return False


##### Quiet mode test
# A few runs (before intermittent mode) were quiet but no valve info
def quiet_mode_check(data):
    quiet_before_valves = ['20251211_5','20251211_6','20251211_7','20251216_13','20260109_2',
                           '20260109_5','20260112_2', '20260113_4','20260113_5','20260115_5']
    if "slow_daq" in data and "valves" in data["slow_daq"] and (data["slow_daq"]["valves"] is not None and not np.all(np.asarray(data["slow_daq"]["valves"]) == 0)):
        valves = data["slow_daq"]["valves"]
        # Check if bit index 2 (third bit) is 0 at any point
        quiet_mode = bool(np.any(((valves >> 2) & 1) == 0))

    elif data['event_info']['run_id'] in quiet_before_valves:
        return True
    else:
        return False
    
    return quiet_mode



def freq_rms(freqs, psd, f_low, f_high):
    mask = (freqs >= f_low) & (freqs < f_high)
    return np.sqrt(np.sum(psd[mask]))





def acoustic_noise(data):
    results = {}

    # Successful expansion test
    results['succ_expansion'] = bool(single_expansion_check_run(data))

    # Quiet mode test
    results['quiet_mode'] = bool(quiet_mode_check(data))


    # Pressure Setpoint
    results['pset'] = data['event_info']['pset_lo']


    wvfs = data['acoustics']['Waveforms'][0]
    ranges_mV = data['acoustics']['Range']

    sample_rate_hz=data['run_control']['acous']['sample_rate']
    sample_rate_hz = float(sample_rate_hz.replace("MS/s", "").strip()) * 1e6

    pre_trig_len=data['run_control']['acous']['pre_trig_len']
    trig_time=pre_trig_len
    post_trig_len=data['run_control']['acous']['post_trig_len']

    
    baselines = []
    rms_noise_full = []
    rms_noise_0_5kHz = []
    rms_noise_5_10kHz = []
    rms_noise_10_20kHz = []
    rms_noise_20_50kHz = []
    rms_noise_50_200kHz = []
    
    piezo_indices = range(1, 7)  # Skip 1 (trigger) and 8 (not working), piezo 5 died (unplugged?) at some point
    for i in piezo_indices:

        wf_adc = np.array(wvfs[i])
        n = len(wf_adc)
        range_mV = ranges_mV[0, i]
        volts_per_count = (range_mV / 1000) / 2 / 32768
        wf = wf_adc * volts_per_count
        

        time_axis = np.arange(n) / sample_rate_hz

        # Baseline and noise pre-trigger
        pretrig_mask = time_axis < trig_time
        pretrig = wf[pretrig_mask]

        baseline = np.mean(pretrig)
        baselines.append(baseline)

        pretrig_centered = pretrig - baseline
        rms = np.sqrt(np.mean((pretrig_centered)**2))
        rms_noise_full.append(rms)



        # Freq bands
        N = len(pretrig_centered)
        fft_vals = np.fft.rfft(pretrig_centered)
        freqs = np.fft.rfftfreq(N, d=1/sample_rate_hz)
    
        psd = np.abs(fft_vals)**2 / N**2
    
        rms_noise_0_5kHz.append(freq_rms(freqs, psd, 0, 5e3))
        rms_noise_5_10kHz.append(freq_rms(freqs, psd, 5e3, 10e3))
        rms_noise_10_20kHz.append(freq_rms(freqs, psd, 10e3, 20e3))
        rms_noise_20_50kHz.append(freq_rms(freqs, psd, 20e3, 50e3))
        rms_noise_50_200kHz.append(freq_rms(freqs, psd, 50e3, 200e3))
    

    results['baselines'] = baselines
    results['rms_noise_full'] = rms_noise_full
    results['rms_noise_0_5kHz'] = rms_noise_0_5kHz
    results['rms_noise_5_10kHz'] = rms_noise_5_10kHz
    results['rms_noise_10_20kHz'] = rms_noise_10_20kHz
    results['rms_noise_20_50kHz'] = rms_noise_20_50kHz
    results['rms_noise_50_200kHz'] = rms_noise_50_200kHz
    
    return results

    