# EventAnalysis.py
import numpy as np

keys = ['run_id', 'event_id', 'ev_exit_code', 'ev_livetime', 'cum_livetime', 'pset_lo', 'pset_hi', 'pset_ramp1', 'pset_ramp_down', 'pset_ramp_up', 'start_time', 'end_time', 'trigger_source']

def EventAnalysis(ev):
    out = dict()
    for k in keys:
        out[k] = np.asarray(ev['event_info'][k])
    return out

if __name__ == "__main__":
    from GetEvent import GetEvent

    TEST_RUN = "/exp/e961/data/SBC-25-daqdata/20250611_1/"
    TEST_EVENT = 0
    print(EventAnalysis(GetEvent(TEST_RUN, TEST_EVENT)))
