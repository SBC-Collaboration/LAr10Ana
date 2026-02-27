#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import os
import matplotlib.pyplot as plt
import SBCcode
from SBCcode.DataHandling import ReadBinary as rb

def pltr(ind):
    sipm = pmttraces['traces'][ind,0] * pmttraces['v_scale'][ind,0] + pmttraces['v_offset'][ind,0]
    trigger = pmttraces['traces'][ind,1] * pmttraces['v_scale'][ind,1] + pmttraces['v_offset'][ind,1]
    plt.figure()
    plt.plot(sipm[:-100])
    plt.plot(trigger[:-100]*0.1)
    plt.show()

def bl(sipm):
    # calculate baseline
    baseline = np.mean(sipm[:200])
    return baseline

def area(ind):
    sipm = pmttraces['traces'][ind,0] * pmttraces['v_scale'][ind,0] + pmttraces['v_offset'][ind,0]
    baseline = bl(sipm)
    # calculate area
    min_ind = np.argmin(sipm[:-200])
    area = np.sum(sipm[min_ind-250:min_ind+350] - baseline) * pmttraces['dt'][ind,0]
    #area = np.sum(sipm[1100:1700] - baseline) * pmttraces['dt'][ind,0]
    return area

def get_bl(ind):
    sipm = pmttraces['traces'][ind,0] * pmttraces['v_scale'][ind,0] + pmttraces['v_offset'][ind,0]
    return bl(sipm)

def get_led(ind):
    trigger = pmttraces['traces'][ind,1] * pmttraces['v_scale'][ind,1] + pmttraces['v_offset'][ind,1]
    # if triggered, more than 10 points will be below -0.1
    return np.sum(trigger[:-200] < -0.1) >10

def get_ind(ind, cut = baseline_cut):
    return np.argwhere(np.cumsum(cut)==ind)[0,0] + 1

raw_directory = "/bluearc/storage/SBC-20-data"

led_runs = ["20201007_0", "20201012_1"]
bg_runs = ["20201007_1", "20201012_0"]
cs_runs = ["20201007_2", "20201007_5"]
ba_runs = ["20201007_3"]
na_runs = ["20201007_4"]
hg_runs = ["20201012_2"]

for run in ["20201015_4"]:
    
    evs = np.array(os.listdir(os.path.join(raw_directory, run)))
    dirfilter = [os.path.isdir(os.path.join(raw_directory, run, i)) for i in evs]
    evs = evs[dirfilter]
    
    all_area = np.array([])
    all_trig_area = np.array([])
    
    for ev in evs:
        pmttraces = rb.ReadBlock(os.path.join(raw_directory, run, str(ev), "PMTtraces.bin"))
        
        # 't0_sec', 't0_frac', 'v_offset', 'v_scale', 't0', 'dt', 'lost_samples', 'traces'
#        t = np.array([pmttraces['t0'][ind,0] + i * pmttraces['dt'][ind,0] for i in range(2048)])
#        sipm = pmttraces['traces'][ind,0] * pmttraces['v_scale'][ind,0] + pmttraces['v_offset'][ind,0]
#        trigger = pmttraces['traces'][ind,1] * pmttraces['v_scale'][ind,1] + pmttraces['v_offset'][ind,1]
        #baseline = np.mean(sipm[:100])
        
        all_bl = np.array([get_bl(ind) for ind in range(len(pmttraces['traces']))])
        
        baseline_lower_cut = (np.mean(all_bl)-np.std(all_bl) <= all_bl)
        baseline_upper_cut = (all_bl <= 0)
        baseline_cut = [baseline_lower_cut[i] and baseline_upper_cut[i] for i in range(len(baseline_lower_cut))]
        
        ev_area = np.array([area(ind) for ind in range(len(pmttraces['traces']))])[baseline_cut]
        all_area = np.append(all_area, ev_area)
#        
#        led_cut = np.array([get_led(ind) for ind in range(len(pmttraces['traces']))])
#        led_and_baseline_cut = [baseline_cut[i] and led_cut[i] for i in range(len(baseline_cut))]
#        ev_trig_area = np.array([area(ind) for ind in range(len(pmttraces['traces']))])[led_and_baseline_cut]
#        all_trig_area = np.append(all_trig_area, ev_trig_area)

#plt.figure()
#plt.hist(-all_area, bins=50, histtype="step", range=(-1e-9, 5e-9))
##plt.hist(-all_trig_area, bins=50, histtype="step", range=(-1e-9, 5e-9))
#plt.ylim(2e2, )
#plt.title("area histogram (Background)")
#plt.xlabel("V*second")
#plt.ylabel("counts")
##plt.legend(["All triggers", " pulsed"])
#plt.yscale("log")
#plt.show()


data = rb.ReadBlock("data.bin")
pixel_factor = 1.4e-9
plt.figure()
#plt.hist(-data["bg_area"]/pixel_factor, bins=500, histtype="step", range=(-1,60))
#plt.hist(-data["cs_area"]/pixel_factor, bins=500, histtype="step", range=(-1,60))
#plt.hist(-data["ba_area"]/pixel_factor, bins=500, histtype="step", range=(-1,60))
#plt.hist(-data["na_area"]/pixel_factor, bins=500, histtype="step", range=(-1,60))
#plt.hist(-data["hg_area"]/pixel_factor, bins=500, histtype="step", range=(-1,60))
plt.hist(-all_area/pixel_factor, bins=50, histtype="step", range=(-1, 60))
plt.title("area histogram")
plt.xlabel("pixels")
plt.ylabel("counts")
#plt.ylim(2e4, )
#plt.legend(["Background", "Cs137", "Ba133", "Na22", "Hg203"])
#plt.legend(["All triggers", " pulsed"])
plt.yscale("log")
#plt.axvline(x=0.85, ls="--", lw=1)
#plt.axvline(x=1.7, ls="--", lw=1)
#plt.axvline(x=2.55, ls="--", lw=1)
plt.show()

#SBCcode.DataHandling.WriteBinary.WriteBinaryNtupleFile("data.bin", data)