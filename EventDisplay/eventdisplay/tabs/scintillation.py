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

        # Reload button
        self.reload_fastdaq_scintillation_button = tk.Button(
            self.scintillation_tab_left,
            text='Reload',
            command=self.draw_fastdaq_scintillation
        )
        self.reload_fastdaq_scintillation_button.grid(row=2, column=0, columnspan=2, sticky='WE')

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
        time = np.arange(len(data)) * (1 / self.fastdaq_event['acoustics']['sample_rate'])

        # Plot
        self.scintillation_ax.clear()
        self.scintillation_ax.plot(time, data)
        self.scintillation_ax.set_title(self.scintillation_combobox.get() + " " + str(self.fastdaq_event["event_info"]["run_id"][0]) + " " + str(self.fastdaq_event["event_info"]["event_id"][0]))
        self.scintillation_ax.set_xlabel('[s]')
        self.scintillation_ax.set_ylabel('[V]')
        self.scintillation_ax.relim()
        self.scintillation_ax.autoscale_view()

        # Render
        self.scintillation_canvas.draw_idle()
        self.scintillation_canvas.get_tk_widget().grid(row=0, column=1, sticky='NW')
