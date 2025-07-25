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
        tk.Frame().__init__(self, master)

        # State variable for loading fastDAQ data
        self.load_fastdaq_scintillation_var = tk.BooleanVar(value=False)

        self.trigger_index = 0

        # Build UI
        self.create_scintillation_widgets()
        self.scintillation_canvas_setup()

    def create_scintillation_widgets(self):
        # Main frame
        self.scintillation_tab = tk.Frame(self.notebook)
        self.notebook.add(self.scintillation_tab, text='SiPM')

        # Left and right panels + Frames to hold buttons and variables
        self.scintillation_tab_left = tk.Frame(self.scintillation_tab, bd=5, relief=tk.SUNKEN)
        self.scintillation_tab_left.grid(row=0, column=0, sticky='NW')
        self.scintillation_tab_right = tk.Frame(self.scintillation_tab, bd=5, relief=tk.SUNKEN)
        self.scintillation_tab_right.grid(row=0, column=1, sticky='NW')
        self.info_frame = tk.LabelFrame(self.scintillation_tab_left, text="Pulse Info", padx=5, pady=5)
        self.info_frame.grid(row=6, column=0, columnspan=4, sticky='WE', pady=(10, 0))
        self.trigger_step_frame = tk.Frame(self.scintillation_tab_left)
        self.trigger_step_frame.grid(row=1, column=4, padx=(10, 0), sticky="W")

        # Load fastDAQ checkbutton
        self.load_fastdaq_scintillation_checkbutton = tk.Checkbutton(
            self.scintillation_tab_left,

            text='Load SiPM',

            variable=self.load_fastdaq_scintillation_var,
            command=self.load_fastdaq_scintillation
        )
        self.load_fastdaq_scintillation_checkbutton.grid(row=0, column=0, columnspan=2, sticky='WE')

        # Channel selector
        tk.Label(self.scintillation_tab_left, text='Channel:').grid(row=1, column=0, sticky='WE')
        self.scintillation_combobox = ttk.Combobox(self.scintillation_tab_left, width=12)

        self.scintillation_combobox.bind("<<ComboboxSelected>>", lambda _: self.new_channel())
        self.scintillation_combobox.grid(row=1, column=1, sticky='WE')
    
        # Label showing number of triggers
        self.trigger_count_label = tk.Label(self.scintillation_tab_left, text="Triggers: ?")
        self.trigger_count_label.grid(row=1, column=2, padx=(10, 0), sticky="W")

        # Entry box for trigger selection
        self.trigger_var = tk.StringVar()
        self.trigger_entry = tk.Entry(self.scintillation_tab_left, textvariable=self.trigger_var, width=6)
        self.trigger_entry.grid(row=1, column=3, padx=(2, 10), sticky="W")
        self.trigger_entry.bind("<Return>", self.on_trigger_entry_change)

        # 6 skip buttons in 3x2 grid

        btn_specs = [
            ("+1", 1),   ("-1", -1),
            ("+10", 10), ("-10", -10),
            ("+100", 100), ("-100", -100)
        ]

        for i, (label, step) in enumerate(btn_specs):
            btn = tk.Button(self.trigger_step_frame, text=label, width=4, command=lambda s=step: self.shift_trigger(s))
            btn.grid(row=i // 2, column=i % 2, padx=1, pady=1)

        # Info box showing pulse characteristics
        self.hit_t0_label = tk.Label(self.info_frame, text="hit_t0: N/A")
        self.hit_t0_label.grid(row=0, column=0, sticky='W')

        self.hit_amp_label = tk.Label(self.info_frame, text="hit_amp: N/A")
        self.hit_amp_label.grid(row=1, column=0, sticky='W')

        self.hit_area_label = tk.Label(self.info_frame, text="hit_area: N/A")
        self.hit_area_label.grid(row=2, column=0, sticky='W')

        self.wvf_area_label = tk.Label(self.info_frame, text="wvf_area: N/A")
        self.wvf_area_label.grid(row=3, column=0, sticky='W')

        self.second_pulse_label = tk.Label(self.info_frame, text="Second Pulse: N/A")
        self.second_pulse_label.grid(row=4, column=0, sticky='W')

        self.baseline_label = ttk.Label(self.info_frame, text="Baseline: ")
        self.baseline_label.grid(row=5, column=0) 

        self.rms_label = ttk.Label(self.info_frame, text="RMS: ")
        self.rms_label.grid(row=6, column=0) 

        self.enable_cuts_var = tk.BooleanVar()
        self.enable_cuts_check = ttk.Checkbutton(
            self.info_frame,
            text="Enable Cuts",
            variable=self.enable_cuts_var,
            command=self.toggle_cut_widgets
        )
        self.enable_cuts_check.grid(row=7, column=0, sticky="w", pady=(10, 0))

        self.cut_widgets = {}

        cut_variables = ["hit_amp", "hit_area", "wvf_area", "hit_t0", "second_pulse"]
        for i, var in enumerate(cut_variables):
            if var == "second_pulse":
                yes_var = tk.BooleanVar(value=True)
                no_var = tk.BooleanVar(value=True)
                yes_check = ttk.Checkbutton(self.info_frame, text="Yes", variable=yes_var)
                no_check = ttk.Checkbutton(self.info_frame, text="No", variable=no_var)

                self.cut_widgets[var] = {
                    "yes_var": yes_var,
                    "no_var": no_var,
                    "yes_check": yes_check,
                    "no_check": no_check
                }
            else:
                entry = ttk.Entry(self.info_frame, width=6)
                mode = tk.StringVar(value="above")
                above_button = ttk.Radiobutton(self.info_frame, text="Above", variable=mode, value="above")
                below_button = ttk.Radiobutton(self.info_frame, text="Below", variable=mode, value="below")

                self.cut_widgets[var] = {
                    "entry": entry,
                    "mode": mode,
                    "above": above_button,
                    "below": below_button
                }

        # Time window slider
        self.t_start_var = tk.DoubleVar(value=0.0)
        self.t_end_var   = tk.DoubleVar(value=0.0)

        tk.Label(self.scintillation_tab_left, text='T start:').grid(row=2, column=0, sticky='E')
        self.t_start_slider = tk.Scale(
            self.scintillation_tab_left,
            variable=self.t_start_var,
            orient='horizontal',
            from_=0, to=0,         
            resolution=1,             
            showvalue=True,
            length=200
        )
        self.t_start_slider.grid(row=2, column=1, sticky='WE')
        self.t_start_slider.bind("<ButtonRelease-1>", self.on_slider_release)


        tk.Label(self.scintillation_tab_left, text='T end:').grid(row=2, column=2, sticky='E')
        self.t_end_slider = tk.Scale(
            self.scintillation_tab_left,
            variable=self.t_end_var,
            orient='horizontal',
            from_=0, to=0,            
            resolution=1,            
            showvalue=True,
            length=200
        )
        self.t_end_slider.grid(row=2, column=3, sticky='WE')
        self.t_end_slider.bind("<ButtonRelease-1>", self.on_slider_release)


            # Voltage sliders
        self.v_lower_var = tk.DoubleVar(value=0.0)
        self.v_upper_var = tk.DoubleVar(value=0.0)

        tk.Label(self.scintillation_tab_left, text='V lower:').grid(row=3, column=0, sticky='E')
        self.v_lower_slider = tk.Scale(
            self.scintillation_tab_left,
            variable=self.v_lower_var,
            orient='horizontal',
            from_=0, to=0,         
            resolution=1e-3,        
            showvalue=True,
            length=200
        )
        self.v_lower_slider.grid(row=3, column=1, sticky='WE')
        self.v_lower_slider.bind("<ButtonRelease-1>", self.on_slider_release)


        tk.Label(self.scintillation_tab_left, text='V upper:').grid(row=3, column=2, sticky='E')
        self.v_upper_slider = tk.Scale(
            self.scintillation_tab_left,
            variable=self.v_upper_var,
            orient='horizontal',
            from_=0, to=0,         
            resolution=1e-3,    
            showvalue=True,
            length=200
        )
        self.v_upper_slider.grid(row=3, column=3, sticky='WE')
        self.v_upper_slider.bind("<ButtonRelease-1>", self.on_slider_release)

        # Frequency cutoff sliders
        self.f_low_var = tk.DoubleVar(value=0.0)
        self.f_high_var = tk.DoubleVar(value=0.0)

        tk.Label(self.scintillation_tab_left, text='F low (Hz):').grid(row=4, column=0, sticky='E')
        self.f_low_slider = tk.Scale(
            self.scintillation_tab_left,
            variable=self.f_low_var,
            orient='horizontal',
            from_=0, to=0,
            resolution=1,
            showvalue=True,
            length=200
        )
        self.f_low_slider.grid(row=4, column=1, sticky='WE')
        self.f_low_slider.bind("<ButtonRelease-1>", self.on_slider_release)

        tk.Label(self.scintillation_tab_left, text='F high (Hz):').grid(row=4, column=2, sticky='E')
        self.f_high_slider = tk.Scale(
            self.scintillation_tab_left,
            variable=self.f_high_var,
            orient='horizontal',
            from_=0, to=0,
            resolution=1,
            showvalue=True,
            length=200
        )
        self.f_high_slider.grid(row=4, column=3, sticky='WE')
        self.f_high_slider.bind("<ButtonRelease-1>", self.on_slider_release)

        # Reload button
        self.reload_fastdaq_scintillation_button = tk.Button(
            self.scintillation_tab_left,
            text='Reload',
            command=self.draw_fastdaq_scintillation
        )
        self.reload_fastdaq_scintillation_button.grid(row=5, column=0, columnspan=4, sticky='WE')


    def scintillation_canvas_setup(self):
        # Figure and canvas for plotting
        self.scintillation_fig = Figure(figsize=(7, 9), dpi=100)
        gs = self.scintillation_fig.add_gridspec(3, 1)

        self.scintillation_ax = self.scintillation_fig.add_subplot(gs[0, 0])
        self.gain_ax = self.scintillation_fig.add_subplot(gs[2, 0])
        self.fft_ax = self.scintillation_fig.add_subplot(gs[1, 0])

        # Spacing and fitting
        self.scintillation_fig.subplots_adjust(left=0.12, right=0.98, top=0.95, bottom=0.05, hspace=0.8)

        self.scintillation_canvas = FigureCanvasTkAgg(self.scintillation_fig, self.scintillation_tab_right)
        self.scintillation_canvas.get_tk_widget().grid(row=0, column=1, sticky='NSEW')

        self.scintillation_tab_right.grid_rowconfigure(0, weight=1)
        self.scintillation_tab_right.grid_columnconfigure(1, weight=1)


    def load_fastdaq_scintillation(self):
        if not self.load_fastdaq_scintillation_var.get():

            self.scintillation_tab_right.grid_forget()
            return
        # Show panel
        self.scintillation_tab_right.grid(row=0, column=1, sticky='NW')

        # Load event data
        selected = ["run_control", "scintillation", "event_info"]
        self.path = os.path.join(self.raw_directory, self.run)
        try:
            self.scint_fastdaq_event = GetEvent(self.path, self.event, *selected)
            self.pulses = SiPMPulses(self.scint_fastdaq_event)
            self.gain = SiPMGain(self.pulses)
            self.photon = PhotonT0(self.pulses)
            # Populate channels
            n_channels = self.scint_fastdaq_event['scintillation']['Waveforms'].shape[1]
            self.scintillation_combobox['values'] = [f"Channel {i+1}" for i in range(n_channels)]
            self.scintillation_combobox.current(0)
            # Initial draw
            self.new_channel()
        except Exception as e:
            print(e)
            self.scint_error()
        # Clean up memory
        gc.collect()

    

    def draw_fastdaq_scintillation(self, val=None):
        if self.scint_fastdaq_event == None:
            self.scint_error()
            return
        if not self.load_fastdaq_scintillation_var.get():
            self.scintillation_tab_right.grid_forget()
            return
        if self.path != os.path.join(self.raw_directory, self.run):
            self.load_fastdaq_scintillation()
        self.scintillation_tab_right.grid(row=0, column=1, sticky='NW')

        # Get range and domain
        start = self.t_start_var.get()
        end   = self.t_end_var.get()
        vlow  = self.v_lower_var.get()
        vhigh = self.v_upper_var.get()
        # Get frequency cutoffs
        flow = self.f_low_var.get()
        fhigh = self.f_high_var.get()
        # Get FFT arrays
        filtered_data = self.filter_signal_by_freq(self.data, flow, fhigh)
        fft_vals = np.fft.rfft(filtered_data)
        freqs = np.fft.rfftfreq(len(filtered_data), d=self.dt)
        fft_mag = np.abs(fft_vals)

        # Plot raw data and filtered data on same axis
        # Plot vertical dotted lines at hit_t0 and hit_end (estimated)
        self.scintillation_ax.clear()
        self.scintillation_ax.plot(self.time, self.data, label='Raw')
        self.scintillation_ax.plot(self.time, filtered_data, label='Filtered')
        channel_index = self.scintillation_combobox.current()
        if channel_index >= 0:
            t0_val = self.pulses['hit_t0'][channel_index, self.trigger_index]
            amp_val = self.pulses['hit_amp'][channel_index, self.trigger_index]
            area_val = self.pulses['hit_area'][channel_index, self.trigger_index]

            if not np.isnan(t0_val) and amp_val > 0:
                self.scintillation_ax.axvline(t0_val, color='red', linestyle='dotted', label='hit_t0')
        # Horizontal line for threshold = baseline + 5 * RMS
        try:
            baseline_val = self.pulses['baseline'][channel_index, self.trigger_index]
            rms_val = self.pulses['rms'][channel_index, self.trigger_index]
            threshold_val = baseline_val + 5 * rms_val
            self.scintillation_ax.axhline(threshold_val, color='purple', linestyle='dotted', label='Threshold')

        except Exception as e:
            print("Failed to add threshold line:", e)
        self.scintillation_ax.relim()
        self.scintillation_ax.autoscale_view()
        self.scintillation_ax.set_xlim(start, end)
        self.scintillation_ax.set_ylim(vlow, vhigh)
        self.scintillation_ax.set_title(self.scintillation_combobox.get() + " Run: " + str(self.run) + " Event: " + str(self.event))
        self.scintillation_ax.set_xlabel('[s]')
        self.scintillation_ax.set_ylabel('[mV]')
        self.scintillation_ax.legend(loc='upper right', fontsize='small')

        # Histogram of hit amplitudes
        amps = self.photon['amp']
        self.gain_ax.clear()
        self.gain_ax.hist(amps[~np.isnan(amps)], bins=50)
        self.gain_ax.set_title("Hits per Amplitude histogram")
        self.gain_ax.set_xlim(0, )
        self.gain_ax.set_xlabel("Pulse amplitude (mV)")
        self.gain_ax.set_ylabel("Hits")

        # Plot FFT
        self.fft_ax.clear()
        self.fft_ax.plot(freqs[1:], fft_mag[1:])  
        self.fft_ax.set_xlim(flow, fhigh)
        self.fft_ax.set_title("FFT Magnitude")
        self.fft_ax.set_xlabel("Frequency (Hz)")
        self.fft_ax.set_ylabel("Magnitude")

        # Update info labels
        channel_index = self.scintillation_combobox.current()
        try:
            hit_t0_val = self.pulses['hit_t0'][channel_index, self.trigger_index]
            hit_amp_val = self.pulses['hit_amp'][channel_index, self.trigger_index]
            hit_area_val = self.pulses['hit_area'][channel_index, self.trigger_index]
            wvf_area_val = self.pulses['wvf_area'][channel_index, self.trigger_index]
            second_pulse_val = self.pulses['second_pulse'][channel_index, self.trigger_index]
            baseline_val = self.pulses['baseline'][channel_index, self.trigger_index]
            rms_val = self.pulses['rms'][channel_index, self.trigger_index]

            self.hit_t0_label.config(text=f"hit_t0: {hit_t0_val:.4e}")
            self.hit_amp_label.config(text=f"hit_amp: {hit_amp_val:.4e}")
            self.hit_area_label.config(text=f"hit_area: {hit_area_val:.4e}")
            self.wvf_area_label.config(text=f"wvf_area: {wvf_area_val:.4e}")
            self.second_pulse_label.config(text=f"Second Pulse: {'Yes' if second_pulse_val else 'No'}")
            self.baseline_label.config(text=f"Baseline: {baseline_val:.4f}")
            self.rms_label.config(text=f"RMS: {rms_val:.4f}")
        except Exception as e:
            print("Error updating pulse info:", e)
            self.hit_t0_label.config(text="hit_t0: N/A")
            self.hit_amp_label.config(text="hit_amp: N/A")
            self.hit_area_label.config(text="hit_area: N/A")
            self.wvf_area_label.config(text="wvf_area: N/A")
            self.second_pulse_label.config(text="Second Pulse: N/A")
            self.baseline_label.config(text="Baseline: N/A")


        # Render
        self.scintillation_canvas.draw_idle()
        self.scintillation_canvas.get_tk_widget().grid(row=0, column=1, sticky='NW')

    def new_channel(self):
        idx = self.scintillation_combobox.current()
        if idx < 0:
            idx = 0
        # Update voltage and time slider range
        self.data = self.scint_fastdaq_event['scintillation']['Waveforms'][self.trigger_index][idx]
        self.time = np.arange(len(self.data)) * (1 / self.scint_fastdaq_event['scintillation']['sample_rate'])
        num_trigs = self.scint_fastdaq_event['scintillation']['Waveforms'].shape[0]
        self.trigger_count_label.config(text=f"Triggers: {num_trigs}")
        self.dt   = self.time[1] - self.time[0] 
        self.t0   = self.time[0]    
        self.t1 = self.time[-1]
        self.v0   = np.min(self.data) 
        self.v1 = np.max(self.data)
        self.dv   = (self.v1 - self.v0) / 100.0 
        self.t_start_slider.config(from_=self.t0, to=self.t1, resolution=self.dt)
        self.t_end_slider.config(from_=self.t0, to=self.t1, resolution=self.dt)
        self.v_lower_slider.config(from_=self.v0, to=self.v1, resolution=self.dv)
        self.v_upper_slider.config(from_=self.v0, to=self.v1, resolution=self.dv)
        self.t_start_slider.set(self.t0)
        self.t_end_slider.set(self.t1)
        self.v_upper_slider.set(self.v1)
        self.v_lower_slider.set(self.v0)

        # Update FFT slider range
        nyquist = 1 / (2 * self.dt)
        self.f_low_slider.config(from_=0, to=nyquist, resolution=nyquist / 100)
        self.f_high_slider.config(from_=0, to=nyquist, resolution=nyquist / 100)
        self.f_low_slider.set(0)
        self.f_high_slider.set(nyquist)

        self.draw_fastdaq_scintillation()

    def scint_error(self):
        self.scint_fastdaq_event = None
        self.scintillation_combobox['values'] = []
        self.scintillation_combobox.set('')
        self.scintillation_ax.clear()
        self.scintillation_ax.text(
            0.5, 0.5,
            f"No Data For Run {self.run} Event {self.event}",
            transform=self.scintillation_ax.transAxes,
            fontsize=10,
            ha='center',
            va='center'
        )
        self.gain_ax.clear()
        self.gain_ax.text(
            0.5, 0.5,
            f"No Gain Data",
            transform=self.gain_ax.transAxes,
            fontsize=10,
            ha='center',
            va='center'
        )
        self.fft_ax.clear()
        self.fft_ax.text(
            0.5, 0.5,
            f"No FFT Data",
            transform=self.fft_ax.transAxes,
            fontsize=10,
            ha='center',
            va='center'
        )
        self.scintillation_canvas.draw_idle()

    def filter_signal_by_freq(self, data, flow, fhigh):
        fft = np.fft.rfft(data)
        freqs = np.fft.rfftfreq(len(data), d=self.dt)

        # Zero out frequencies outside the desired band
        bandpass = (freqs >= flow) & (freqs <= fhigh)
        fft[~bandpass] = 0

        # Inverse FFT to return filtered time signal
        return np.fft.irfft(fft, n=len(data))

    # Wrapper function for sliders
    def on_slider_release(self, var):
        self.draw_fastdaq_scintillation()   

    def on_trigger_entry_change(self, event):
        try:
            idx = int(self.trigger_var.get())
            max_idx = self.scint_fastdaq_event['scintillation']['Waveforms'].shape[0]
            if 0 <= idx < max_idx:
                self.trigger_index = idx
                self.new_channel()
            else:
                print(f"Trigger index {idx} out of range.")
        except ValueError:
            print("Invalid trigger index entered.")

    def shift_trigger(self, step):
        max_idx = self.scint_fastdaq_event['scintillation']['Waveforms'].shape[0]
        new_idx = self.trigger_index + step
        new_idx = max(0, min(new_idx, max_idx - 1))  # Clamp to valid range

        self.trigger_index = new_idx
        self.trigger_var.set(str(new_idx))
        self.new_channel()

    def toggle_cut_widgets(self):
        if self.enable_cuts_var.get():
            for i, (var, widgets) in enumerate(self.cut_widgets.items()):
                row = 0 + i
                if var == "second_pulse":
                    widgets["yes_check"].grid(row=row, column=3)
                    widgets["no_check"].grid(row=row, column=4)
                else:
                    widgets["entry"].grid(row=row, column=3)
                    widgets["above"].grid(row=row, column=4)
                    widgets["below"].grid(row=row, column=5)
        else:
            for var, widgets in self.cut_widgets.items():
                if var == "second_pulse":
                    widgets["yes_check"].grid_remove()
                    widgets["no_check"].grid_remove()
                else:
                    widgets["entry"].grid_remove()
                    widgets["above"].grid_remove()
                    widgets["below"].grid_remove()
