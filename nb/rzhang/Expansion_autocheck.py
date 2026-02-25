from sbcbinaryformat import Streamer, Writer
import numpy as np
import matplotlib.pyplot as plt

from GetEvent import GetEvent

from ana import AcousticT0 
from scipy.signal import firwin, filtfilt
from scipy.optimize import least_squares
import importlib

importlib.reload(AcousticT0)


def single_expansion_check_run(path,i):
    TEST_RUN = "/exp/e961/data/SBC-25-daqdata/20250611_1/"
    TEST_RUN2 = "/exp/e961/app/users/runze/data/20251120_12/"
    # TEST_RUN2 = "/exp/e961/data/users/gputnam/SBC-25-daqdata-test/20251103_1/"
    TEST_EVT = i
    
    data = GetEvent(path, i,strictMode=False)
    data.keys()
    
    data["cam"].keys()
    
    data["acoustics"].keys()
    
    data["acoustics"]["Range"]

    pre_trig_len=data['run_control']['acous']['pre_trig_len']
    print("acoustic pretrigger post trigger in s",pre_trig_len,"\n",pre_trig_len,data['run_control']['acous']['post_trig_len'])
    pre_trig_frames=data['run_control']['cams']['cam1']['buffer_len']-data['run_control']['cams']['cam1']['post_trig']
    post_trig_frames=data['run_control']['cams']['cam1']['post_trig']
    print("camera pretrigger/post trigger frames",pre_trig_frames, post_trig_frames)
    
    
    wvfs = data["acoustics"]["Waveforms"]
    wvfs_psi = wvfs*(-35/2**16)

    piezo0 = wvfs_psi[0, 7, :]

    average_window = 100
    start_pressure = np.mean(piezo0[:average_window])
    end_pressure = np.mean(piezo0[600000:600000+average_window])
    if start_pressure-end_pressure<0.3:
        return True
    else:
        return False



def single_expansion_data_run(path,i):
    # this return the difference of pressure at 0 and 600ms
    TEST_RUN = "/exp/e961/data/SBC-25-daqdata/20250611_1/"
    TEST_RUN2 = "/exp/e961/app/users/runze/data/20251120_12/"
    # TEST_RUN2 = "/exp/e961/data/users/gputnam/SBC-25-daqdata-test/20251103_1/"
    TEST_EVT = i
    
    data = GetEvent(path, i,strictMode=False)
    data.keys()
    
    data["cam"].keys()
    
    data["acoustics"].keys()
    
    data["acoustics"]["Range"]

    pre_trig_len=data['run_control']['acous']['pre_trig_len']
    print("acoustic pretrigger post trigger in s",pre_trig_len,"\n",pre_trig_len,data['run_control']['acous']['post_trig_len'])
    pre_trig_frames=data['run_control']['cams']['cam1']['buffer_len']-data['run_control']['cams']['cam1']['post_trig']
    post_trig_frames=data['run_control']['cams']['cam1']['post_trig']
    print("camera pretrigger/post trigger frames",pre_trig_frames, post_trig_frames)
    
    
    wvfs = data["acoustics"]["Waveforms"]
    wvfs_psi = wvfs*(-35/2**16)

    piezo0 = wvfs_psi[0, 7, :]

    average_window = 100
    start_pressure = np.mean(piezo0[:average_window])
    end_pressure = np.mean(piezo0[600000:600000+average_window])
    return start_pressure-end_pressure
        

# parameter_list is your event list mask. True means success expansion in the event, False means failure.
parameter_list = []
for i in range(5,6,1):
# You can change i range in case you want to loop in diffrent data set
    try:
        base_path = f"/exp/e961/app/users/runze/data/20251120_{i}/"
        for j in range(0,99,1):
            # j is the event number in one dataset
            try:
                print(i,j, "new looop", base_path)
                paras = single_expansion_check_run(base_path,j)
                parameter_list.append(paras)
            except:
                continue
    except:
        continue
print(parameter_list)




































































































