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

def single_expansion_slowdaq_check_run(path,i):
    TEST_RUN = "/exp/e961/data/SBC-25-daqdata/20250611_1/"
    TEST_RUN2 = "/exp/e961/app/users/runze/data/20251120_12/"
    # TEST_RUN2 = "/exp/e961/data/users/gputnam/SBC-25-daqdata-test/20251103_1/"
    TEST_EVT = i
    
    data = GetEvent(path, i,strictMode=False)
    data.keys()
    
    compress_idx = 0
    for i in range(len(data["slow_daq"]['SERVO3321_OUT'])):
        if data["slow_daq"]['SERVO3321_OUT'][i]>75:
            print(i,data["slow_daq"]['SERVO3321_OUT'][i],data["slow_daq"]['time_ms'][i])
            compress_time =data["slow_daq"]['time_ms'][i]
            compress_idx = i
            break

            
   
    print("compress",compress_time)

    # find t_start
    time_cut  = -10
    slowdaq_PT = data["slow_daq"]['PT2121']
    slowdaq_time = data["slow_daq"]['time_ms']
    t_start_time_window = slowdaq_time[:compress_idx+time_cut]
    pressure_cut_compress = slowdaq_PT[:compress_idx+time_cut]
    reverse_time_start = t_start_time_window[::-1]
    reverse_time_start = [int(x) for x in reverse_time_start]
    reverse_pressure =pressure_cut_compress[::-1]
    time_width = reverse_time_start[-2]-reverse_time_start[-1]
    p_set = float(data["event_info"]['pset_hi'])
    # print("set", p_set, "mode",p_mode)
    if reverse_pressure[0]<p_set+0.05:
        for i in range(len(reverse_pressure)):
            if reverse_pressure[i]<p_set+0.05 and reverse_pressure[i+1]>p_set+0.05:
                interpolation_part = (p_set+0.05-reverse_pressure[i])*(reverse_time_start[i+1]-reverse_time_start[i])/(reverse_pressure[i+1]-reverse_pressure[i])
                interpolation_time = reverse_time_start[i] + interpolation_part
                # print((p_set+0.05-reverse_pressure[i]),(reverse_time_start[i+1]-reverse_time_start[i]),(reverse_pressure[i+1]-reverse_pressure[i]))
                # print(interpolation_part, interpolation_time)
                expansion_time = interpolation_time
        
        
                break
       
                
    else:
        expansion_time = compress_time
    

    print("expansion",expansion_time) 
    time_diff = compress_time-expansion_time
    if time_diff>1000: # 1st for compression and expansion
        return (True, time_diff)
    else:
        return (False,time_diff)

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




































































































