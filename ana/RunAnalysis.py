# RunAnalysis.py
import numpy as np

keys = ['run_id', 'run_exit_code', 'num_events', 'run_livetime', 'comment', 'run_start_time', 'run_end_time', 'active_modules', 'pset_mode', 'pset_lo', 'pset_hi', 'source1_ID', 'source1_location', 'rc_ver', 'red_caen_ver', 'niusb_ver', 'sbc_binary_ver']

def RunAnalysis(ev):
    out = dict()
    run_info = ev["run_info"]

    # convert the run info into new format
    if ("pset_lo" not in run_info or "pset_hi" not in run_info) and "pset" in run_info:
        run_info["pset_lo"] = run_info["pset"]
        run_info["pset_hi"] = run_info["pset"]
    
    for k in keys:
        value = run_info.get(k)
        out[k] = np.asarray(value) if value is not None else np.array([])
    return out

if __name__ == "__main__":
    from GetEvent import GetEvent

    TEST_RUN = "/exp/e961/data/SBC-25-daqdata/20250611_1/"
    TEST_EVENT = 0
    print(RunAnalysis(GetEvent(TEST_RUN, TEST_EVENT)))
