import numpy as np

SAMPLE_FREQ = 62.5 # MHz

def SiPMPulses(ev):
    # Stuff to save, with defaults
    default_output = dict(
        baseline=np.array([]),
        rms=np.array([]),
        hit_t0=np.array([]),
        hit_area=np.array([]),
        hit_amp=np.array([]),
        wvf_area=np.array([]),
        second_pulse=np.array([]),
        wvf_timestamp=np.array([]),
    )
    out = default_output

    if ev is None:
        return out

    # One event has N SiPMs readout M times, each with T samples
    # Each readout is stored is ev["sipm_traces"], with a shape (M, N, T)

    # For each SiPM, for each readout, obtain the t0 and the voltage

    # Waveform in (ticks, ADC)
    traces = ev["scintillation"]["Waveforms"].T.astype(float)

    # Subtract offset and convert to mV
    for i_sipm in range(traces.shape[1]):
        group = i_sipm // 8
        chan  = i_sipm % 8
        group_ctrl = ev["run_control"]["caen"]["group%i" % group]
        offset = group_ctrl["offset"] + group_ctrl["ch_offset"][chan]
        range_mV = float(group_ctrl["range"][:-4])

        traces[:, i_sipm, :] -= offset
        traces[:, i_sipm, :] *= range_mV


    decimation = ev["run_control"]["caen"]["global"]["decimation"]
    sample_rate =  SAMPLE_FREQ/(2**decimation)

    # obtain the leading baseline and RMS
    N_SAMPLE_BASELINE = 40

    baseline = traces[:N_SAMPLE_BASELINE].mean(axis=0)
    rms = traces[:N_SAMPLE_BASELINE].std(axis=0)

    # flip the trace and correct for baseline
    trace_V = -(traces - baseline)

    # Start time of hit
    N_SIGMA_THRESHOLD = 5
    above_threshold = trace_V > rms*N_SIGMA_THRESHOLD
    t0_ind = np.argmax(above_threshold, axis=0)
    t0 = t0_ind / sample_rate

    # Final time of hit
    tf_ind = np.argmax(np.cumsum(~above_threshold, axis=0) > t0_ind, axis=0)

    # build index into each waveform
    wvf_index = np.zeros(trace_V.shape, dtype=np.int32)
    wvf_index[:, :, :] = np.arange(0, wvf_index.shape[0]).reshape((wvf_index.shape[0], 1, 1))

    # mask area inside hit
    hit_trace_V = trace_V*(wvf_index >= t0_ind)*(wvf_index < tf_ind)

    # voltage values
    hit_area = hit_trace_V.sum(axis=0)
    hit_amplitude = hit_trace_V.max(axis=0)
    wvf_area = (trace_V*(wvf_index >= t0_ind)).sum(axis=0)

    # Are there any secondary hits after the first one?
    N_SECONDPULSE_TICK_DELAY = 10
    second_pulse = (above_threshold*(wvf_index > (tf_ind + N_SECONDPULSE_TICK_DELAY))).any(axis=0).astype(int)

    out["baseline"] = baseline
    out["rms"] = rms
    out["hit_t0"] = t0
    out["hit_area"] = hit_area
    out["hit_amp"] = hit_amplitude
    out["wvf_area"] = wvf_area
    out["second_pulse"] = second_pulse

    # If no hit was found, set relevant values to nan
    out["hit_t0"][t0 == 0] = np.nan
    out["hit_area"][t0 == 0] = np.nan
    out["hit_amp"][t0 == 0] = np.nan
    out["wvf_area"][t0 == 0] = np.nan
    out["second_pulse"][t0 == 0] = False

    return out

if __name__ == "__main__":
    from GetEvent import GetEvent

    TEST_RUN = "/exp/e961/data/SBC-25-daqdata/20250611_1/"
    TEST_EVENT = 0
    out = SiPMPulses(GetEvent(TEST_RUN, TEST_EVENT))
    print("OUTPUTS:", *out.keys())
    print("BASELINES:", ", ".join(map(str, list(out["baseline"].mean(axis=-1)))))
    print("RMS:", ", ".join(map(str, list(out["rms"].mean(axis=-1)))))
    print("NHIT:", ", ".join(map(str, list((~np.isnan(out["hit_t0"])).astype(int).sum(axis=-1)))))
