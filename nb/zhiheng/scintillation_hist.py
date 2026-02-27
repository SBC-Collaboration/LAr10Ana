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
    endpoint = 4e-7
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


alpha_runs = ["20201106_0", "20201109_1"]
bg_runs = ["20201007_1", "20201012_0"]
led_runs = ["20201111_0"]
cs_runs = ["20201007_2", "20201007_5"]
cs_runs_new = ["20201113_0"]
alpha_runs_low_threshold = ["20201113_1"]
alpha_extended = ["20201116_1"]
cs_extended = ["20201015_4", "20201016_5"]
cs_extended_new = ["20201116_2"]
bg_extended = ["20201016_4"]
alpha_high_gain = ["20201116_3"]

bg55 = ["20201117_1"]
bg51 = ["20201117_3"]
bg53 = ["20201117_5"]
bg49 = ["20201117_7"]
cs55 = ["20201117_2"]
cs51 = ["20201117_4"]
cs53 = ["20201117_6"]
cs49 = ["20201117_8"]

gas_bg = ["20201123_0"]
gas_cs = ["20201123_1"]
gas_bg_55 = ["20201123_2"]
gas_bg_53 = ["20201123_3"]
gas_bg_51 = ["20201123_4"]
gas_cs_55 = ["20201123_5"]
gas_cs_53 = ["20201123_6"]
gas_cs_51 = ["20201123_7"]
vacuum_bg_55 = ["20201123_8"]
vacuum_bg_53 = ["20201123_9"]
vacuum_bg_51 = ["20201123_10"]

no_alpha_gas_55 = ["20201203_0", "20201203_3"]
no_alpha_gas_53 = ["20201203_1"]
no_alpha_gas_51 = ["20201203_2"]
liquid_bg55 = ["20201205_0"]
liquid_bg53 = ["20201205_1"]
liquid_bg51 = ["20201205_2"]
liquid_cs55 = ["20201205_3"]
liquid_cs53 = ["20201205_4"]
liquid_cs51 = ["20201205_5"]

#bg_area, bg_livetime = gen_hist(bg_runs)
#alpha_area, alpha_livetime = gen_hist(alpha_runs)
#led_area, led_livetime = gen_hist(led_runs)
#cs_area, cs_livetime = gen_hist(cs_runs)
#cs_area_new, cs_livetime_new = gen_hist(cs_runs_new)
#alpha_low_area, alpha_low_livetime = gen_hist(alpha_runs_low_threshold)
alpha_ext_area, alpha_ext_livetime = gen_hist(alpha_extended)
cs_ext_new_area, cs_ext_new_livetime = gen_hist(cs_extended_new)
cs_ext_area, cs_ext_livetime = gen_hist(cs_extended)
bg_ext_area, bg_ext_livetime = gen_hist(bg_extended)
#alpha_high_gain_area, alpha_high_gain_livetime = gen_hist(alpha_high_gain)
#bg55_area, bg55_livetime = gen_hist(bg55)
#bg53_area, bg53_livetime = gen_hist(bg53)
#bg51_area, bg51_livetime = gen_hist(bg51)
#bg49_area, bg49_livetime = gen_hist(bg49)
#cs55_area, cs55_livetime = gen_hist(cs55)
#cs53_area, cs53_livetime = gen_hist(cs53)
#cs51_area, cs51_livetime = gen_hist(cs51)
#cs49_area, cs49_livetime = gen_hist(cs49)
#gas_bg_area, gas_bg_livetime = gen_hist(gas_bg)
#gas_cs_area, gas_cs_livetime = gen_hist(gas_cs)
gas_bg55_area, gas_bg55_livetime = gen_hist(gas_bg_55)
gas_bg53_area, gas_bg53_livetime = gen_hist(gas_bg_53)
gas_bg51_area, gas_bg51_livetime = gen_hist(gas_bg_51)
gas_cs55_area, gas_cs55_livetime = gen_hist(gas_cs_55)
#gas_cs53_area, gas_cs53_livetime = gen_hist(gas_cs_53)
#gas_cs51_area, gas_cs51_livetime = gen_hist(gas_cs_51)
vacuum_bg55_area, vacuum_bg55_livetime = gen_hist(vacuum_bg_55)
#vacuum_bg53_area, vacuum_bg53_livetime = gen_hist(vacuum_bg_53)
#vacuum_bg51_area, vacuum_bg51_livetime = gen_hist(vacuum_bg_51)
no_alpha_55_area, no_alpha_55_livetime = gen_hist(no_alpha_gas_55)
#no_alpha_53_area, no_alpha_53_livetime = gen_hist(no_alpha_gas_53)
#no_alpha_51_area, no_alpha_51_livetime = gen_hist(no_alpha_gas_51)
liquid_bg55_area, liquid_bg55_livetime = gen_hist(liquid_bg55)
liquid_bg53_area, liquid_bg53_livetime = gen_hist(liquid_bg53)
liquid_bg51_area, liquid_bg51_livetime = gen_hist(liquid_bg51)
liquid_cs55_area, liquid_cs55_livetime = gen_hist(liquid_cs55)

#plot_hist(bg_area, bg_livetime)
#plot_hist(alpha_area, alpha_livetime)
##plot_hist(led_area, 180)
#plot_hist(cs_area, cs_livetime)
#plot_hist(cs_area_new, cs_livetime_new)
#plot_hist(alpha_low_area, alpha_low_livetime)
#plt.legend(["background", "alpha", "Cs", "Cs with alpha", "alpha (low threshold)", "diff in Cs", "diff in background"])

#extended plots
#plt.figure()
#plot_hist(bg_ext_area, bg_ext_livetime)
#plot_hist(alpha_ext_area, alpha_ext_livetime)
#plot_hist(cs_ext_area, cs_ext_livetime)
#plot_hist(cs_ext_new_area, cs_ext_new_livetime)
#plt.legend(["background", "alpha", "Cs", "Cs with alpha"])

#extended plots
plt.figure()
plot_hist(liquid_bg55_area, liquid_bg55_livetime)
plot_hist(bg_ext_area, bg_ext_livetime)
plot_hist(alpha_ext_area, alpha_ext_livetime)
plt.legend(["Old Bg", "New Bg", "Alpha"])
#plot_hist(cs_ext_area, cs_ext_livetime)
#plot_hist(liquid_cs55_area, liquid_cs55_livetime)
#plot_hist(cs_ext_new_area, cs_ext_new_livetime)
#plt.legend(["Old Cs", "New Cs", "Alpha+Cs"])
#plt.title("Cs Optical Geometry Comparison")

endpoint = 4e-7
binwidth = (endpoint+1e-9)/100
p = np.histogram(-no_alpha_55_area, bins=100, range=(-1e-9, endpoint))
x = p[1][:-1]
y1 = p[0]/no_alpha_55_livetime/binwidth
plt.step(x,y1)
p = np.histogram(-gas_bg55_area, bins=100, range=(-1e-9, endpoint))
x = p[1][:-1]
y2 = p[0]/gas_bg55_livetime/binwidth
plt.step(x,y2)
diff1 = y2-y1
plt.step(x, diff1)
plt.xlabel("V*second")
plt.ylabel("rate (Hz/V/s)")
#plt.ylim(15, 3e5)
plt.yscale("log")
plt.legend(["Bg", "Bg+Alpha", "Diff"])
plt.title("Gas Background Difference")

p = np.histogram(-liquid_bg55_area, bins=100, range=(-1e-9, endpoint))
x = p[1][:-1]
y1 = p[0]/liquid_bg55_livetime/binwidth
plt.step(x,y1)
p = np.histogram(-alpha_ext_area, bins=100, range=(-1e-9, endpoint))
x = p[1][:-1]
y2 = p[0]/alpha_ext_livetime/binwidth
plt.step(x,y2)
diff2 = y2-y1
plt.step(x, diff)
plt.xlabel("V*second")
plt.ylabel("rate (Hz/V/s)")
#plt.ylim(15, 3e5)
plt.yscale("log")
plt.legend(["Bg", "Bg+Alpha", "Diff"])
plt.title("Liquid Background Difference")

endpoint = 4e-7
binwidth = (endpoint+1e-9)/100
p = np.histogram(-liquid_cs55_area, bins=100, range=(-1e-9, endpoint))
x = p[1][:-1]
y1 = p[0]/liquid_cs55_livetime/binwidth
plt.step(x,y1)
p = np.histogram(-cs_ext_new_area, bins=100, range=(-1e-9, endpoint))
x = p[1][:-1]
y2 = p[0]/cs_ext_new_livetime/binwidth
plt.step(x,y2)
diff3 = y2-y1
plt.step(x, diff3)
plt.xlabel("V*second")
plt.ylabel("rate (Hz/V/s)")
#plt.ylim(15, 3e5)
plt.yscale("log")
plt.legend(["Cs", "Cs+Alpha", "Diff"])
plt.title("Gas Background Difference")

plt.step(np.flip(x), np.cumsum(np.flip(diff1))*binwidth)
plt.step(np.flip(x), np.cumsum(np.flip(diff2))*binwidth)
plt.step(np.flip(x), np.cumsum(np.flip(np.abs(diff3)))*binwidth)
plt.legend(["Bg-liquid", "Bg-gas", "Cs-liquid"])
plt.xlabel("V*second")
plt.ylabel("rate (Hz)")
plt.title("Reverse CDF of Alpha rates with Background subtracted")