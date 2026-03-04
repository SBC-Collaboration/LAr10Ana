# RunAnalysis.py
import numpy as np

keys = ['run_id', 'run_exit_code', 'num_events', 'run_livetime', 'comment', 'run_start_time', 'run_end_time', 'active_modules', 'pset_mode', 'pset_lo', 'pset_hi', 'source1_ID', 'source1_location', 'rc_ver', 'red_caen_ver', 'niusb_ver', 'sbc_binary_ver']

def RunAnalysis(ev):
    out = dict()
    for k in keys:
        out[k] = np.asarray(ev['run_info'][k])
    return out

if __name__ == "__main__":
    from GetEvent import GetEvent

    TEST_RUN = "/exp/e961/data/SBC-25-daqdata/20250611_1/"
    TEST_EVENT = 0
    print(RunAnalysis(GetEvent(TEST_RUN, TEST_EVENT)))
