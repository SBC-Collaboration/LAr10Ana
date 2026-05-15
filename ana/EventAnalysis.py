# EventAnalysis.py
import numpy as np

keys = ['run_id', 'event_id', 'ev_exit_code', 'ev_livetime', 'cum_livetime', 'pset_lo', 'pset_hi', 'pset_ramp1', 'pset_ramp_down', 'pset_ramp_up', 'start_time', 'end_time', 'trigger_source']

def EventAnalysis(ev):
    out = dict()
    event_info = ev["event_info"]

    # convert the run info into new format
    if ("pset_lo" not in event_info or "pset_hi" not in event_info) and "pset" in event_info:
        event_info["pset_lo"] = event_info["pset"]
        event_info["pset_hi"] = event_info["pset"]
    
    for k in keys:
        value = event_info.get(k)
        out[k] = np.asarray(value) if value is not None else np.array([])
    return out

if __name__ == "__main__":
    from GetEvent import GetEvent

    TEST_RUN = "/exp/e961/data/SBC-25-daqdata/20250611_1/"
    TEST_EVENT = 0
    print(EventAnalysis(GetEvent(TEST_RUN, TEST_EVENT)))
