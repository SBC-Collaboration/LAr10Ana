# Gary's scintillation rate and # of hit SiPMs analysis code

# Import Libraries
import sys
import numpy as np

# Functions to run code

# FFT cuts based on FFT shape
# Good waveforms should have max FFT near 0 MHz
def _FFT_frequency_filtering(waveforms, dt, freq_cutoff_hz):
    # Number of samples in waveform
    N = waveforms.shape[-1]

    # Compute the FFT values and frequencies
    fft_vals = np.abs(np.fft.fft(waveforms, axis=-1))
    fft_freq = np.abs(np.fft.fftfreq(N, d=dt))

    # Find where the maximum FFT val occurs
    max_amp_indices = np.argmax(fft_vals[:, :, 1:], axis=-1)

    # Don't need the fft_vals anymore
    del fft_vals

    # Get the frequency a which the max FFT val occurs
    peak_frequencies_hz = fft_freq[max_amp_indices]

    # Create a mask for noisy waveforms
    bad_channel_mask = peak_frequencies_hz >= freq_cutoff_hz

    # If it's a noisy waveform, replace it with zeros.
    filtered_waveforms = waveforms.copy()
    filtered_waveforms[bad_channel_mask] = 0.
    return filtered_waveforms


# Computing the signal strength of waveforms
def _calculate_signal_strength(fft_filtered_waveforms):

    # Compute the baseline of each waveform
    baseline = np.mean(fft_filtered_waveforms, axis=2)

    # Calculate signal strength
    signalPeak = np.min(fft_filtered_waveforms, axis=2) - baseline
    signalDroop = np.max(fft_filtered_waveforms, axis=2) - baseline
    signalStrength = signalDroop + signalPeak

    return signalStrength

# Simple function to unwrap CAEN timestamps
def _unwrap_caen_timestamp(ts, max_ts):
    ts = np.asarray(ts, dtype=np.int64)

    # Detect rollovers
    rollovers = np.diff(ts, axis=-1, prepend=0) < 0

    # Cumulative count of rollovers
    rollover_count = np.cumsum(rollovers, axis=-1)
    return ts + rollover_count * max_ts

# Main function to compute SiPMs hit per each file
def ScintillationRateAnalysis(ev):
    # default value
    output = {
        "n_hits": np.zeros(1, dtype=np.uint8),
        "good_event_mask": np.zeros(1, dtype=np.uint32),
    }
    if ev is None or not ev['event_info']['loaded'] or not ev['scintillation']['loaded']:
        print("File not loaded. Quitting.")
        return output

    scint = ev["scintillation"]
    # Load the waveforms
    waveforms = scint['Waveforms']()

    # load other data which may be important
    sample_rate = 62.5e6
    decimation = ev['run_control']['caen']['global']['decimation']
    scint_timestamps = _unwrap_caen_timestamp(scint['TriggerTimeTag'](), 2**31)
    livetime = scint_timestamps[-1] * 8e-9  # timestmap is 8 ns

    if decimation == 0: decimation = 1
    dt = 1 / (sample_rate * decimation)        

    # Denoise waveforms and collect signal strength
    fft_filtered_waveforms = _FFT_frequency_filtering(
        waveforms=waveforms, 
        dt=dt, 
        freq_cutoff_hz = 0.6 # MHz (Still rough guess, works well)
    )
    
    signal_strength = _calculate_signal_strength(fft_filtered_waveforms)

    # Find all signals with strength < -10 (should be good pulses).
    signal_strength_limit = -10.
    mask = signal_strength < signal_strength_limit

    # Sum the total number of hits in the event
    NHits = np.sum(mask, axis=1).astype(np.uint8)
    output["n_hits"] = NHits

    # convert boolean mask array to uint32
    weights = 2**np.arange(32, dtype=np.uint32)
    output["good_event_mask"] = mask.dot(weights)

    return output

# The rest of this code can be used in a jupyter notebook to see the results. 
if __name__ == "__main__":
    # from GetEvent import GetEvent
    ana_path = "../LAr10Ana/"
    sys.path.insert(0, ana_path)
    import GetEvent

    # Compute the scintillation rates
    test_run = "/exp/e961/data/SBC-25-daqdata/20251107_28.tar"
    test_event = 0
    background_ev =  GetEvent.GetEvent(test_run, test_event, "run_control", "run_info", "event_info", "scintillation", 
                                strictMode=False, lazy_load_scintillation=False)
    out = ScintillationRateAnalysis(background_ev)
    for k in out.keys():
        print(f"{k}: {out[k][:5]}")
