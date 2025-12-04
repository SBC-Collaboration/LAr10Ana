import numpy as np

keys = ['run_id', 'event_id', 'ev_exit_code', 'ev_livetime', 'cum_livetime', 'pset', 'pset_hi', 'pset_slope', 'pset_period', 'start_time', 'end_time', 'trigger_source']

def EventAnalysis(ev):
    default_output = dict()
    for k in keys:
        default_output[k] = np.asarray(-1.)

    try:
        out = default_output
        for k in keys:
            out[k] = np.asarray(ev['event_info'][k])
        return out
    except:
        return default_output

if __name__ == "__main__":
    from GetEvent import GetEvent

    TEST_RUN = "/exp/e961/data/SBC-25-daqdata/20250611_1/"
    TEST_EVENT = 0
    print(EventAnalysis(GetEvent(TEST_RUN, TEST_EVENT)))
