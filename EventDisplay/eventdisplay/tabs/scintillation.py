import gc
import os
import sys
import numpy as np
import matplotlib
import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Hacky
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
from GetEvent import GetEvent
# Even more hacky
BASE = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                     '..', '..', '..'))
ANA_DIR = os.path.join(BASE, 'ana')
sys.path.insert(0, ANA_DIR)
from SiPMPulses import SiPMPulses
from SiPMGain import SiPMGain
from PhotonT0 import PhotonT0

class Scintillation(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)

        # State variable for loading fastDAQ data
        self.load_fastdaq_scintillation_var = tk.BooleanVar(value=False)

        # Build UI
        self.create_scintillation_widgets()
        self.scintillation_canvas_setup()

    def create_scintillation_widgets(self):
        # Main frame
        self.scintillation_tab = tk.Frame(self.notebook)
        self.notebook.add(self.scintillation_tab, text='SiPM')

        # Left and right panels
        self.scintillation_tab_left = tk.Frame(self.scintillation_tab, bd=5, relief=tk.SUNKEN)
        self.scintillation_tab_left.grid(row=0, column=0, sticky='NW')
        self.scintillation_tab_right = tk.Frame(self.scintillation_tab, bd=5, relief=tk.SUNKEN)
        self.scintillation_tab_right.grid(row=0, column=1, sticky='NW')

        # Load fastDAQ checkbutton
        self.load_fastdaq_scintillation_checkbutton = tk.Checkbutton(
            self.scintillation_tab_left,
            text='Load fastDAQ',
            variable=self.load_fastdaq_scintillation_var,
            command=self.load_fastdaq_scintillation
        )
        self.load_fastdaq_scintillation_checkbutton.grid(row=0, column=0, columnspan=2, sticky='WE')

        # Channel selector
        tk.Label(self.scintillation_tab_left, text='Channel:').grid(row=1, column=0, sticky='WE')
        self.scintillation_combobox = ttk.Combobox(self.scintillation_tab_left, width=12)
        self.scintillation_combobox.grid(row=1, column=1, sticky='WE')
    
        # Time window slider
        self.t_start_var = tk.DoubleVar(value=0.0)
        self.t_end_var   = tk.DoubleVar(value=0.0)

        tk.Label(self.scintillation_tab_left, text='T start:').grid(row=2, column=0, sticky='E')
        self.t_start_slider = tk.Scale(
            self.scintillation_tab_left,
            variable=self.t_start_var,
            orient='horizontal',
            from_=0, to=0,            # will be updated on reload
            resolution=1,             # will be updated on reload
            showvalue=True,
            length=200
        )
        self.t_start_slider.grid(row=2, column=1, sticky='WE')

        tk.Label(self.scintillation_tab_left, text='T end:').grid(row=2, column=2, sticky='E')
        self.t_end_slider = tk.Scale(
            self.scintillation_tab_left,
            variable=self.t_end_var,
            orient='horizontal',
            from_=0, to=0,            # will be updated on reload
            resolution=1,             # will be updated on reload
            showvalue=True,
            length=200
        )
        self.t_end_slider.grid(row=2, column=3, sticky='WE')

            # ── Voltage‐window sliders ──
        self.v_lower_var = tk.DoubleVar(value=0.0)
        self.v_upper_var = tk.DoubleVar(value=0.0)

        # V lower
        tk.Label(self.scintillation_tab_left, text='V lower:').grid(row=3, column=0, sticky='E')
        self.v_lower_slider = tk.Scale(
            self.scintillation_tab_left,
            variable=self.v_lower_var,
            orient='horizontal',
            from_=0, to=0,            # will update on reload
            resolution=1e-3,          # will update on reload
            showvalue=True,
            length=200
        )
        self.v_lower_slider.grid(row=3, column=1, sticky='WE')

        # V upper
        tk.Label(self.scintillation_tab_left, text='V upper:').grid(row=3, column=2, sticky='E')
        self.v_upper_slider = tk.Scale(
            self.scintillation_tab_left,
            variable=self.v_upper_var,
            orient='horizontal',
            from_=0, to=0,            # will update on reload
            resolution=1e-3,          # will update on reload
            showvalue=True,
            length=200
        )
        self.v_upper_slider.grid(row=3, column=3, sticky='WE')



        # Reload button
        self.reload_fastdaq_scintillation_button = tk.Button(
            self.scintillation_tab_left,
            text='Reload',
            command=self.draw_fastdaq_scintillation
        )
        self.reload_fastdaq_scintillation_button.grid(row=4, column=0, columnspan=4, sticky='WE')

    def scintillation_canvas_setup(self):
        # Figure and canvas for plotting
        self.scintillation_fig = Figure(figsize=(7, 5), dpi=100)
        self.scintillation_ax = self.scintillation_fig.add_subplot(111)
        self.scintillation_canvas = FigureCanvasTkAgg(self.scintillation_fig, self.scintillation_tab_right)

    def load_fastdaq_scintillation(self):
        if not self.load_fastdaq_scintillation_var.get():
            # Hide when unchecked
            self.scintillation_tab_right.grid_forget()
            return
        # Show panel
        self.scintillation_tab_right.grid(row=0, column=1, sticky='NW')

        # Load event data
        selected = ["run_control", "scintillation", "event_info"]
        self.path = os.path.join(self.raw_directory, self.run)
        self.fastdaq_event = GetEvent(self.path, self.event, *selected)
        self.pulses = SiPMPulses(self.fastdaq_event)
        self.gain = SiPMGain(self.pulses)
        self.photon = PhotonT0(self.pulses)


        # Populate channels
        n_channels = self.fastdaq_event['scintillation']['Waveforms'].shape[1]
        self.scintillation_combobox['values'] = [f"Channel {i+1}" for i in range(n_channels)]

        # Initial draw
        self.draw_fastdaq_scintillation()

        # Clean up memory
        gc.collect()

    def draw_fastdaq_scintillation(self):
        if not self.load_fastdaq_scintillation_var.get():
            self.scintillation_tab_right.grid_forget()
            return
        self.scintillation_tab_right.grid(row=0, column=1, sticky='NW')

        # Retrieve selected channel data
        idx = self.scintillation_combobox.current()
        if idx < 0:
            return
        data = self.fastdaq_event['scintillation']['Waveforms'][0][idx]
        time = np.arange(len(data)) * (1 / self.fastdaq_event['scintillation']['sample_rate'])

        # … after loading data and time …
        dt   = time[1] - time[0] if len(time) > 1 else 1.0
        t0   = time[0];    t1 = time[-1]
        v0   = np.min(data); v1 = np.max(data)
        dv   = (v1 - v0) / 100.0 if v1 != v0 else 1e-3

        # update time sliders (you already have this)
        self.t_start_slider.config(from_=t0, to=t1, resolution=dt)
        self.t_end_slider  .config(from_=t0, to=t1, resolution=dt)

        # update voltage sliders
        self.v_lower_slider.config(from_=v0, to=v1, resolution=dv)
        self.v_upper_slider.config(from_=v0, to=v1, resolution=dv)

        # update slider ranges & step
        self.t_start_slider.config(from_=t0, to=t1, resolution=dt)
        self.t_end_slider.config(from_=t0, to=t1, resolution=dt)

        # Plot
        self.scintillation_ax.clear()
        self.scintillation_ax.plot(time, data)
        self.scintillation_ax.relim()
        self.scintillation_ax.autoscale_view()

        start = self.t_start_var.get();    end   = self.t_end_var.get()
        vlow  = self.v_lower_var.get();    vhigh = self.v_upper_var.get()

        self.scintillation_ax.set_xlim(start, end)
        self.scintillation_ax.set_ylim(vlow, vhigh)
        self.scintillation_ax.set_title(self.scintillation_combobox.get() + " Run: " + str(self.fastdaq_event["event_info"]["run_id"][0]) + " Event: " + str(self.fastdaq_event["event_info"]["event_id"][0]))
        self.scintillation_ax.set_xlabel('[s]')
        self.scintillation_ax.set_ylabel('[V]')

        # Render
        self.scintillation_canvas.draw_idle()
        self.scintillation_canvas.get_tk_widget().grid(row=0, column=1, sticky='NW')
