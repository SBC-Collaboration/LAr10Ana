import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import warnings


# Creating function for histograms
def histogram(channels, datatype, title, lb=0, rb=0):
    leftbound = lb
    rightbound = rb

    # Check that bounds are valid
    if not (leftbound == 0 and rightbound == 0) and leftbound >= rightbound:
        print(f"Invalid bounds: leftbound ({leftbound}) >= rightbound ({rightbound})")
        return  # Exit early if bounds are incorrect

    plt.figure(figsize=(8, 5))

    for chan in channels:
        plt.hist(datatype[chan, :], bins=300, histtype='step', label=f"Channel {chan}")

    if not (leftbound == 0 and rightbound == 0):
        plt.axvline(leftbound, c='r', alpha=0.4, label="Lower Bound")
        plt.axvline(rightbound, c='r', alpha=0.4, label="Upper Bound")

    plt.xlabel(f'mV')
    plt.ylabel("Count")
    plt.title(f'{title} Histogram')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()




# Creating function for plotting
def visualcheck(dataType, ind, chan, title = '', fillAll=False):
    plt.figure(figsize=(10, 5))
    
    raw_wf = dataType['raw_wvfm'][ind, chan, :]
    smooth_wf = dataType['smooth_wvfm'][ind, chan, :]

    plt.plot(raw_wf, label="Raw Waveform", alpha=0.7)
    plt.plot(smooth_wf, label="Smoothed Waveform", linewidth=2)

    plt.axvline(dataType['hit_t0_ind'][chan, ind], c='r')
    plt.axvline(dataType['hit_tf_ind'][chan, ind], c='r')

    if title == '':
        plt.title(f"Channel {chan} - Waveform {ind}")
    else: 
        plt.title(f"{title} - Channel {chan} - Waveform {ind}")
    plt.xlabel("Sample")
    plt.ylabel("mV")
    if not fillAll:
        plt.ylim(-50, 150)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()