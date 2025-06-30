import tkinter as tk
from pandas import DataFrame
import gc
import os
import numpy as np
import time
import matplotlib.pyplot as plt
from time import sleep
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from tkinter import ttk, messagebox, DISABLED, NORMAL
from matplotlib.figure import Figure
from tkinter import *
import sys
from PIL import Image, ImageTk
import scipy.signal
import matplotlib
import matplotlib.pyplot as plt
# matplotlib.use('TkAgg')
matplotlib.use('Agg')
from scipy.signal import butter, sosfilt
import getpass
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
from GetEvent import GetEvent

class Scintillation(tk.Frame):
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        
        self.scintillation_cutoff_low = 2000
        self.scintillation_cutoff_high = 10000
        self.scintillation_beginning_time = -.1
        self.scintillation_ending_time = 0.0
        self.incremented_scintillation_event = False
        self.scintillation_timerange_checkbutton_var = tk.BooleanVar(value=False)

        self.create_scintillation_widgets()
        self.scintillation_canvas_setup()

    def scintillation_canvas_setup(self):
        self.scintillation_fig = Figure(figsize=(7, 5), dpi=100)
        self.scintillation_ax = self.scintillation_fig.add_subplot(111)
        self.scintillation_canvas = FigureCanvasTkAgg(self.scintillation_fig, self.scintillation_tab_right)
    
    def create_scintillation_widgets(self):
        self.scintillation_tab = tk.Frame(self.notebook)
        self.notebook.add(self.scintillation_tab, text='scintillation')

        # scintillations tab
        # First setup frames for scintillations tab
        self.scintillation_tab_left = tk.Frame(self.scintillation_tab, bd=5, relief=tk.SUNKEN)
        self.scintillation_tab_left.grid(row=0, column=0, sticky='NW')

        self.scintillation_tab_right = tk.Frame(self.scintillation_tab, bd=5, relief=tk.SUNKEN)
        self.scintillation_tab_right.grid(row=0, column=1, sticky='NW')

        #         self.scintillation_scrollbar = tk.Scrollbar(self.scintillation_tab_right, orient = 'vertical')
        #         self.scintillation_scrollbar.pack(side = 'left', fill = 'y')
        #         self.scintillation_scrollbar.grid(row = 0, column = 0, sticky = tk.N + tk.S + tk.W + tk.E)

        # Now within the scintillations frames setup stuff
        self.load_fastDAQ_scintillation_checkbutton = tk.Checkbutton(
            self.scintillation_tab_left,
            text='Load fastDAQ',
            variable=self.load_fastDAQ_piezo_checkbutton_var,
            command=self.load_fastDAQ_scintillation)
        self.load_fastDAQ_scintillation_checkbutton.grid(row=0, column=0, columnspan=2, sticky='WE')

        self.scintillation_label = tk.Label(self.scintillation_tab_left, text='scintillation:')
        self.scintillation_label.grid(row=1, column=0, sticky='WE')

        self.scintillation_combobox = ttk.Combobox(self.scintillation_tab_left, width=12)
        self.scintillation_combobox.grid(row=1, column=1, sticky='WE')

        self.scintillation_cutoff_low_label = tk.Label(self.scintillation_tab_left, text='Freq cutoff low:')
        self.scintillation_cutoff_low_label.grid(row=2, column=0, sticky='WE')

        self.scintillation_cutoff_low_entry = tk.Entry(self.scintillation_tab_left, width=12)
        self.scintillation_cutoff_low_entry.insert(0, self.scintillation_cutoff_low)
        self.scintillation_cutoff_low_entry.grid(row=2, column=1, sticky='WE')

        self.scintillation_cutoff_high_label = tk.Label(self.scintillation_tab_left, text='Freq cutoff high:')
        self.scintillation_cutoff_high_label.grid(row=3, column=0, sticky='WE')

        self.scintillation_cutoff_high_entry = tk.Entry(self.scintillation_tab_left, width=12)
        self.scintillation_cutoff_high_entry.insert(0, self.scintillation_cutoff_high)
        self.scintillation_cutoff_high_entry.grid(row=3, column=1, sticky='WE')

        self.scintillation_timerange_checkbutton = tk.Checkbutton(
            self.scintillation_tab_left, text='Full time window',
            variable=self.scintillation_timerange_checkbutton_var,
            command=self.draw_fastDAQ_scintillation)
        self.scintillation_timerange_checkbutton.grid(row=6, column=0, columnspan=2, sticky='WE')

        self.scintillation_beginning_time_label = tk.Label(self.scintillation_tab_left, text='Beginning Time:')
        self.scintillation_beginning_time_label.grid(row=4, column=0, sticky='WE')

        self.scintillation_beginning_time_entry = tk.Entry(self.scintillation_tab_left, width=12)
        self.scintillation_beginning_time_entry.insert(0, self.scintillation_beginning_time)
        self.scintillation_beginning_time_entry.grid(row=4, column=1, sticky='WE')

        self.scintillation_ending_time_label = tk.Label(self.scintillation_tab_left, text='Ending Time:')
        self.scintillation_ending_time_label.grid(row=5, column=0, sticky='WE')

        self.scintillation_ending_time_entry = tk.Entry(self.scintillation_tab_left, width=12)
        self.scintillation_ending_time_entry.insert(0, self.scintillation_ending_time)
        self.scintillation_ending_time_entry.grid(row=5, column=1, sticky='WE')

        self.scintillation_plot_t0_checkbutton = tk.Checkbutton(
            self.scintillation_tab_left,
            text='Show t0',
            variable=self.piezo_plot_t0_checkbutton_var,
            command=self.draw_fastDAQ_scintillation)
        self.scintillation_plot_t0_checkbutton.grid(row=7, column=0, columnspan=2, sticky='WE')

        self.reload_fastDAQ_scintillation_button = tk.Button(self.scintillation_tab_left, text='reload',
                                                     command=self.draw_fastDAQ_scintillation)
        self.reload_fastDAQ_scintillation_button.grid(row=8, column=0, sticky='WE')

    def draw_fastDAQ_scintillation(self):
        return
    
    def load_fastDAQ_scintillation(self):
        return