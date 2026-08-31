# Gary's scintillation rate and # of hit SiPMs analysis code

# Import Libraries
import sys
import numpy as np
from ana.BatchSiPMs import BatchSiPMs

# Functions to run code

# Unwrap the CAEN timestamp to account for rollovers
def _unwrap_caen_timestamp(ts, max_ts):
    ts = np.asarray(ts, dtype=np.int64)

    # Detect rollovers
    rollovers = np.diff(ts) < 0

    # Cumulative count of rollovers
    rollover_count = np.concatenate([[0], np.cumsum(rollovers)])
    return ts + rollover_count * max_ts

# Calculate trigger rate using 1s window around this trigger
def _trigger_rate_windowed(ev, window_sec=1.0):
    timetag = ev["scintillation"]['TriggerTimeTag']
    if callable(timetag):
        timetag = timetag()
    timestamp = _unwrap_caen_timestamp(timetag, 2**31) * 8.e-9
    half_win = window_sec / 2.0
    left_idx = np.searchsorted(timestamp, timestamp - half_win, side='left')
    right_idx = np.searchsorted(timestamp, timestamp + half_win, side='right')
    rate = (right_idx - left_idx) / window_sec
    return rate[:, np.newaxis]

# Calculate the ratio of the SiPM pulse to the baseline
def _signal_ratio_filtering(waveforms):

    # Compute the baseline of each waveform
    baseline = np.mean(waveforms[:,:,:50], axis=-1)
    waveforms = waveforms - baseline[:, :, np.newaxis]

    # Calculate the ratio of the signal, noise should be close to 1 (cancels out)
    signalDown = np.min(waveforms, axis=-1)
    signalUp = np.max(waveforms, axis=-1)
    with np.errstate(divide='ignore', invalid='ignore'):
        signalRatio = np.abs(signalDown/signalUp)
    signalRatio = np.nan_to_num(signalRatio, nan=0.0, posinf=0.0, neginf=0.0)

    return signalRatio

def ScintillationRateBatched(ev, nwvf_batch=1000, maxwvf=-1, progress=False, njob=1):
    output = BatchSiPMs(ev, ScintillationRateAnalysis, nwvf_batch=nwvf_batch, maxwvf=maxwvf, progress=progress, njob=njob)
    if ev is not None and ev['scintillation']['loaded']:
        output["trigger_rate"] = _trigger_rate_windowed(ev)[:output["n_hits"].shape[0]].astype(np.uint16)

    return output

# BatchSipm for hit counting only. Trigger rate calculation is done separately in the above function
def ScintillationRateAnalysis(ev):
    # default value
    output = {
        "n_hits": np.zeros((1, 1), dtype=np.uint8),
        "hits_mask": np.zeros((1, 1), dtype=np.uint32),
        "trigger_rate": np.zeros((1,1), dtype=np.uint16)
    }
    if ev is None:
        # Used by BatchSipm to get data structure
        return output
    elif not ev['scintillation']['loaded']:
        print("File not loaded. Quitting.")
        return output
    # Load the waveforms
    waveforms = ev["scintillation"]['Waveforms']
    # Remove non-functional SiPMs
    nonfunctional_sipms = [24, 31]
    waveforms[:, nonfunctional_sipms, :] = 0
    
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
    ana_path = "../LAr10Ana/"
    sys.path.insert(0, ana_path)
    import GetEvent

    # Compute the scintillation rates
    test_run = "/exp/e961/data/SBC-25-daqdata/20251107_28.tar"
    test_event = 0
    background_ev = GetEvent.GetEvent(test_run, test_event, "run_control", "run_info", "event_info", "scintillation",
                                       strictMode=False, lazy_load_scintillation=False)
    # Use ScintillationRateBatched (not ScintillationRateAnalysis directly) to also get the trigger_rate key.
    out = ScintillationRateBatched(background_ev)
    for k in out.keys():
        print(f"{k}: {out[k][:5]}")
