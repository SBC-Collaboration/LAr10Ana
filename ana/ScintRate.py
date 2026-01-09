# Gary's scintillation rate and # of hit SiPMs analysis code

# Import Libraries
import sys
import numpy as np

# Functions to run code

# Calculate the ratio of the SiPM pulse to the baseline
def _signal_ratio_filtering(waveforms):

    # Compute the baseline of each waveform
    baseline = np.mean(waveforms[:,:,:50], axis=-1)
    waveforms = waveforms - baseline[:, :, np.newaxis]

    # Calculate the ratio of the signal, noise should be close to 1 (cancels out)
    signalDown = np.min(waveforms, axis=-1)
    signalUp = np.max(waveforms, axis=-1)
    signalRatio = np.abs(signalDown/signalUp)

    return signalRatio

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
        "n_hits": np.zeros((1,1), dtype=np.uint8),
        "hits_mask": np.zeros((1,1), dtype=np.uint32),
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
    
    signal_ratio = _signal_ratio_filtering(waveforms)

    # Find all signals with a ratio of peak to baseline > 3, these should be good pulses.
    signal_ratio_limit = 3.
    mask = signal_ratio >= signal_ratio_limit

    # Sum the total number of hits in the event
    NHits = np.sum(mask, axis=1).astype(np.uint8)
    output["n_hits"] = NHits[:, np.newaxis]

    # convert boolean mask array to uint32
    weights = 2**np.arange(32, dtype=np.uint32)
    output["hits_mask"] = mask.dot(weights)[:, np.newaxis]

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
