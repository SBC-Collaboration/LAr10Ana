import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import SBCcode
from SBCcode.DataHandling import ReadBinary as rb

def area(ind, pmttraces):
    sipm = pmttraces['traces'][ind,0] * pmttraces['v_scale'][ind,0] + pmttraces['v_offset'][ind,0]
    baseline = np.mean(sipm[:200])
    # calculate area
    min_ind = np.argmin(sipm[:-200])
    area = np.sum(sipm[min_ind-250:min_ind+350] - baseline) * pmttraces['dt'][ind,0]
    #area = np.sum(sipm[1100:1700] - baseline) * pmttraces['dt'][ind,0]
    return area

def get_bl(ind, pmttraces):
    sipm = pmttraces['traces'][ind,0] * pmttraces['v_scale'][ind,0] + pmttraces['v_offset'][ind,0]
    return np.mean(sipm[:200])

def get_led(ind, pmttraces):
    trigger = pmttraces['traces'][ind,1] * pmttraces['v_scale'][ind,1] + pmttraces['v_offset'][ind,1]
    # if triggered, more than 10 points will be below -0.1
    return np.sum(trigger[:-200] < -0.1) >10

raw_directory = "/bluearc/storage/SBC-20-data"

def gen_hist_led(runs):
    all_area = np.array([])
    all_trig_area = np.array([])
    livetime = 0
    for run in runs:
        evs = np.array(os.listdir(os.path.join(raw_directory, run)))
        dirfilter = [os.path.isdir(os.path.join(raw_directory, run, i)) for i in evs]
        evs = evs[dirfilter]
        
        for ev in evs:
            temp = pd.read_csv(os.path.join(raw_directory, "20201109_0", str(0), "Event.txt"), header=None, sep=" ")
            livetime += np.array(temp)[0,-1]
            
            pmttraces = rb.ReadBlock(os.path.join(raw_directory, run, str(ev), "PMTtraces.bin"))          
            all_bl = np.array([get_bl(ind, pmttraces) for ind in range(len(pmttraces['traces']))])
            
            baseline_lower_cut = (np.mean(all_bl)-np.std(all_bl) <= all_bl)
            baseline_upper_cut = (all_bl <= 0)
            baseline_cut = [baseline_lower_cut[i] and baseline_upper_cut[i] for i in range(len(baseline_lower_cut))]
            
            ev_area = np.array([area(ind, pmttraces) for ind in range(len(pmttraces['traces']))])[baseline_cut]
            all_area = np.append(all_area, ev_area)
            
            led_cut = np.array([get_led(ind, pmttraces) for ind in range(len(pmttraces['traces']))])
            led_and_baseline_cut = [baseline_cut[i] and led_cut[i] for i in range(len(baseline_cut))]
            ev_trig_area = np.array([area(ind, pmttraces) for ind in range(len(pmttraces['traces']))])[led_and_baseline_cut]
            all_trig_area = np.append(all_trig_area, ev_trig_area)
    
    plt.figure()
    plt.hist(-all_area, bins=50, histtype="step", range=(-1e-9, 5e-9))
    plt.hist(-all_trig_area, bins=50, histtype="step", range=(-1e-9, 5e-9))
    plt.title("area histogram")
    plt.xlabel("V*second")
    plt.ylabel("counts")
    plt.legend(["All triggers", " pulsed"])
    plt.yscale("log")
    plt.show()

def gen_hist(runs):
    all_area = np.array([])
    livetime = 0
    for run in runs:
        evs = np.array(os.listdir(os.path.join(raw_directory, run)))
        dirfilter = [os.path.isdir(os.path.join(raw_directory, run, i)) for i in evs]
        evs = evs[dirfilter]
        
        for ev in evs:            
            pmttraces = rb.ReadBlock(os.path.join(raw_directory, run, str(ev), "PMTtraces.bin"))
            tstart = np.mean(pmttraces['t0_sec'][0]++pmttraces['t0_frac'][0])
            tend = np.mean(pmttraces['t0_sec'][-1]++pmttraces['t0_frac'][-1])
            livetime += tend - tstart
            
            all_bl = np.array([get_bl(ind, pmttraces) for ind in range(len(pmttraces['traces']))])
            
            baseline_lower_cut = (np.mean(all_bl)-np.std(all_bl) <= all_bl)
            baseline_upper_cut = (all_bl <= 0)
            baseline_cut = [baseline_lower_cut[i] and baseline_upper_cut[i] for i in range(len(baseline_lower_cut))]
            
            ev_area = np.array([area(ind, pmttraces) for ind in range(len(pmttraces['traces']))])[baseline_cut]
            all_area = np.append(all_area, ev_area)
    print("done: ", livetime)
    return all_area, livetime
    
def plot_hist(all_area, livetime):
    #plt.figure()
    #p = plt.hist(-all_area, bins=100, histtype="step", range=(-1e-9, 9e-8))
    endpoint = 1e-7
    binwidth = (endpoint+1e-9)/100
    p = np.histogram(-all_area, bins=100, range=(-1e-9, endpoint))
    x = p[1][:-1]
    y = p[0]/livetime/binwidth
    plt.step(x,y)
    plt.xlabel("V*second")
    plt.ylabel("rate (Hz/V/s)")
    #plt.ylim(15, 3e5)
    plt.yscale("log")
    #plt.show()
    
gas_bg = ["20201123_0"]
gas_bg_area, gas_bg_livetime = gen_hist(gas_bg)
plot_hist(gas_bg_area, gas_bg_livetime)