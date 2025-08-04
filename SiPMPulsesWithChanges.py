import numpy as np
SAMPLE_FREQ = 62.5  # MHz

def SiPMPulses(ev, N_avg = 3):
    
    # Configuration
    N_AVG = N_avg               # Number of samples to average in smoothing (must be odd)
    N_SAMPLE_BASELINE = 50  # Number of samples to average for baseline
    N_SIGMA_THRESHOLD = 2 #-np.inf   # Threshold for hit detection
    n = 10                  # indices used for area window

    # Stuff to save, with defaults
    default_output = dict(
        raw_wvfm=np.array([]),
        smooth_wvfm=np.array([]),
        baseline=np.array([]),
        rms=np.array([]),
        hit_t0_ind=np.array([]),
        hit_tf_ind=np.array([]),
        hit_t0=np.array([]),
        hit_tf=np.array([]),
        hit_area=np.array([]),
        hit_amp=np.array([]),
        wvf_area=np.array([]),
        wvf_timestamp=np.array([]),
        bl_health={}
    )
    out = default_output
    if ev is None:
        return out

    # One event has N SiPMs readout M times, each with T samples
    traces = ev["scintillation"]["Waveforms"].T.astype(float)

    # Subtract offset and convert to mV
    for i_sipm in range(traces.shape[1]):
        group = i_sipm // 8
        chan = i_sipm % 8
        group_ctrl = ev["run_control"]["caen"]["group%i" % group]
        offset = group_ctrl["offset"] + group_ctrl["ch_offset"][chan]
        range_mV = float(group_ctrl["range"][:-4])
        traces[:, i_sipm, :] -= offset
        traces[:, i_sipm, :] *= range_mV

    decimation = ev["run_control"]["caen"]["global"]["decimation"]
    sample_rate = SAMPLE_FREQ / (2 ** decimation)

    # Baseline and RMS
    baseline = traces[:N_SAMPLE_BASELINE].mean(axis=0)
    rms = traces[:N_SAMPLE_BASELINE].std(axis=0)

    # Flip and subtract baseline
    trace_V = -(traces - baseline)

    # Smooth waveform using N_AVG-point moving average
    pad = N_AVG // 2
    padded_trace = np.pad(trace_V, ((pad, pad), (0, 0), (0, 0)), mode='edge')
    smooth_trace = np.zeros_like(trace_V)

    for i in range(N_AVG):
        smooth_trace += padded_trace[i:i + trace_V.shape[0]]
    smooth_trace /= N_AVG

    # Hit finding
    t0_ind = np.argmax(smooth_trace, axis=0) - n
    tf_ind = np.argmax(smooth_trace, axis=0) + int(2.5 * n)

    t0 = t0_ind / SAMPLE_FREQ
    tf = tf_ind / SAMPLE_FREQ

    wvf_index = np.zeros(smooth_trace.shape, dtype=np.int32)
    wvf_index[:, :, :] = np.arange(0, wvf_index.shape[0]).reshape((smooth_trace.shape[0], 1, 1))

    hit_trace_V = smooth_trace * (wvf_index >= t0_ind) * (wvf_index < tf_ind)
    hit_area = hit_trace_V.sum(axis=0)
    hit_amplitude = hit_trace_V.max(axis=0)
    wvf_area = (smooth_trace * (wvf_index >= t0_ind)).sum(axis=0)

    # Zero out results below noise threshold
    
    hit_mask = hit_amplitude < (rms * N_SIGMA_THRESHOLD)
    t0_ind[hit_mask] = 0
    tf_ind[hit_mask] = 0
    t0[hit_mask] = np.nan
    tf[hit_mask] = np.nan
    hit_area[hit_mask] = np.nan
    hit_amplitude[hit_mask] = np.nan
    wvf_area[hit_mask] = np.nan

    # Final outputs
    out["raw_wvfm"] = trace_V.transpose(2, 1, 0)
    out["smooth_wvfm"] = smooth_trace.transpose(2, 1, 0)
    out["baseline"] = baseline
    out["rms"] = rms
    out["hit_t0"] = t0
    out["hit_tf"] = tf
    out["hit_area"] = hit_area
    out["hit_amp"] = hit_amplitude
    out["wvf_area"] = wvf_area
    out["hit_t0_ind"] = t0_ind
    out["hit_tf_ind"] = tf_ind



    
    # Compute and attach baseline health dictionary
    out["bl_health"] = compute_bl_health(rms, baseline)

    return out

# ── Helper Function for Baseline Health ───────────────────────
def compute_bl_health(rms, baseline):
    n_chan, n_evt = rms.shape

    bl_health = {
        "rms_mean": np.zeros((n_chan, 1)),
        "rms_std": np.zeros((n_chan, 1)),
        "rms_outliers": [None] * n_chan,

        "bl_mean": np.zeros((n_chan, 1)),
        "bl_std": np.zeros((n_chan, 1)),
        "bl_outliers": [None] * n_chan
    }

    for ch in range(n_chan):
        # Process RMS
        rms_vals = rms[ch, :]
        rms_mean = np.mean(rms_vals)
        rms_std = np.std(rms_vals)
        rms_outlier_indices = np.where(np.abs(rms_vals - rms_mean) > 3 * rms_std)[0]

        bl_health["rms_mean"][ch, 0] = rms_mean
        bl_health["rms_std"][ch, 0] = rms_std
        bl_health["rms_outliers"][ch] = rms_outlier_indices.tolist()


        
        # Process Baseline
        bl_vals = baseline[ch, :]
        bl_mean = np.mean(bl_vals)
        bl_std = np.std(bl_vals)
        bl_outlier_indices = np.where(np.abs(bl_vals - bl_mean) > 3 * bl_std)[0]

        bl_health["bl_mean"][ch, 0] = bl_mean
        bl_health["bl_std"][ch, 0] = bl_std
        bl_health["bl_outliers"][ch] = bl_outlier_indices.tolist()

    return bl_health





# ── Test usage ──────────────────────────────────────────────
if __name__ == "__main__":
    from GetEvent import GetEvent
    TEST_RUN = "/exp/e961/data/SBC-25-daqdata/20250611_1/"
    TEST_EVENT = 0
    out = SiPMPulses(GetEvent(TEST_RUN, TEST_EVENT))

    print("OUTPUTS:", *out.keys())
    print("BASELINES:", ", ".join(map(str, list(out["baseline"].mean(axis=-1)))))
    print("RMS:", ", ".join(map(str, list(out["rms"].mean(axis=-1)))))
    print("NHIT:", ", ".join(map(str, list((~np.isnan(out["hit_t0"])).astype(int).sum(axis=-1)))))
