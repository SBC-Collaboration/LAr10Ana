import numpy as np

def EventAnalysis(ev):
    default_output = dict()
    try:
        if not ev['event_info']['loaded']:
            return default_output
        out = default_output
        for k in ev['event_info']:
            if not (k == 'loaded'):
                out['Event_' + k] = np.asarray(ev['event_info'][k])
                
        return out
    except:
        return default_output

if __name__ == "__main__":
    from GetEvent import GetEvent

    TEST_RUN = "/exp/e961/data/SBC-25-daqdata/20250611_1/"
    TEST_EVENT = 0
    print(EventAnalysis(GetEvent(TEST_RUN, TEST_EVENT)))
