from sbcbinaryformat import Streamer, Writer
import numpy as np
import matplotlib.pyplot as plt

from GetEvent import GetEvent

from ana import AcousticT0
from scipy.signal import firwin, filtfilt
import importlib
importlib.reload(AcousticT0)

class PT_t0_Finding():
    def __init__(self):
        self.path = "/exp/e961/data/users/gputnam/SBC-25-daqdata-test/20251103_1/"
        self.t0_list = []
        self.event_length = 10
        self.plot_path = "/somepath"
    def main_body(self):
        for i in range(self.event_length):
            self.loop_event(i)

    def plot_t0_histgram(self):
        # lin_bins = np.arange(0, 1e5, 50)
        bin_num  =50

        counts_ini, bin_edges_ini, patches_ini = plt.hist(self.t0_list, bins=bin_num)
        plt.xlabel("Energy (eV)", fontsize=16)
        plt.ylabel(r"Rate (event/hr)", fontsize=16)
        plt.savefig(self.plot_path)
    def loop_event(self,i): # loop around event number
        data = GetEvent(self.path, i, strictMode=False)
        self.find_t0(data)


    def find_t0(self, data):
        wvfs = data["acoustics"]["Waveforms"]

        # plot the first triggered waveform in each channel
        # if channel 7 raw reading is ADC -35/1e4
        # if channel 7 raw reading is bits, -35/2**15
        # MB.PT1101 := AIn.PT1101 * (-35.0* el3052bits);
        # wvfs_psi = wvfs*(35/2**15)
        wvfs_psi = wvfs * (35 / 1e4)
        (time_list_ms, piezoslope0, piezo0_filtered) = self.filter_noise(wvfs_psi)
        t0 = self.quadratic_fitting(time_list_ms, piezoslope0, piezo0_filtered)
        self.t0_list.append(t0)

    def filter_noise(self,wave): # filtered electric noise to get clean waveform
        piezo0 = wave[0, 7, :]
        xlimit = [220, 320]
        ylimit = [-22.5, -20]
        # time in miliseconds
        total_time = len(piezo0)
        time_list_ms = [i / 1e3 for i in range(0, total_time, 1)]

        # add low pass filter
        # assuming 1 microsecond time resolution 1e6Hz
        numtaps = 5000  # filter length (longer = sharper cutoff)
        Fs = 1000000  # sampling rate
        Fc = 10  # low pass filter in Hz
        fir = firwin(numtaps, Fc, window='hamming', fs=Fs)

        piezo0_filtered = filtfilt(fir, [1.0], piezo0)
        piezoslope0 = -piezo0_filtered[1:] + piezo0_filtered[:-1]

        return (time_list_ms, piezoslope0, piezo0_filtered)
    def quadratic_fitting(self, time_list_ms, piezoslope0, piezo0_filtered):
        # quadratic fittings
        # finding ending point, hard cut, find first time when the rate > 1e-4 psi/ms
        # starting point, from [10:]
        # then fitting find t0

        starting_indx = 10000
        # starting 10ms after data collection
        print(piezoslope0[starting_indx])
        ending_indx = 0
        hardcut_threshold = 1e-4
        for i in range(starting_indx, len(piezoslope0), 1):
            if piezoslope0[i] > hardcut_threshold:
                ending_indx = i

                break
        # print("index",ending_indx)
        pressure_before_fit = piezo0_filtered[starting_indx:ending_indx]
        time_fit_range = time_list_ms[starting_indx:ending_indx]

        # qua fit
        coeffs = np.polyfit(time_fit_range, pressure_before_fit, 2)
        a, b, c = coeffs
        print(a, b, c)

        # find t0
        t0_fitted = round(-b / (2 * a), 3)

        return t0_fitted


