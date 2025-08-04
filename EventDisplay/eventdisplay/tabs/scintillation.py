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
from GetEvent import GetEvent, GetScint
from GetEvent import GetEvent, GetScint
# Even more hacky
BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
ANA_DIR = os.path.join(BASE, 'ana')
sys.path.insert(0, ANA_DIR)
from SiPMPulses import SiPMPulses, SiPMPulsesBatched
from SiPMPulses import SiPMPulses, SiPMPulsesBatched
from PhotonT0 import PhotonT0

class Scintillation(tk.Frame):
    def __init__(self, master=None):
        tk.Frame().__init__(self, master)

        # State variable for loading fastDAQ data
        self.load_fastdaq_scintillation_var = tk.BooleanVar(value=False)

        self.trigger_index = 0

        self.channel_info_labels = {}

        self.analysis_cache = None
        self.analyzed_triggers = None


        self.analysis_cache = None
        self.analyzed_triggers = None


        # Build UI
        self.create_scintillation_widgets()
        self.scintillation_canvas_setup()

    # Loading files 
    def load_fastdaq_scintillation(self):
        # Return if the checkbox is not selected
        if not self.load_fastdaq_scintillation_var.get():
            self.scintillation_tab_right.grid_forget()
            return
        self.scintillation_tab_right.grid(row=0, column=1, sticky='NW')
        self.path = os.path.join(self.raw_directory, self.run)
        try:
            self.load_event()
            self.populate_channel_listbox()
            self.reset_analysis_cache()
            self.new_channel()
        # Error handling (prints out error but not where error is. Uncomment to debug)
        except Exception as e:
            print(e)
            self.scint_error()
        gc.collect()

    # Load bin and sbc files
    def load_event(self):
        selected = ["run_control", "scintillation", "event_info"]
        self.scint_fastdaq_event = GetEvent(self.path, self.event, *selected, lazy_load_scintillation=False)

    # Create channel names in listbox
    def populate_channel_listbox(self):
        wf_full = self.scint_fastdaq_event["scintillation"]["Waveforms"]
        self.scintillation_listbox.delete(0, tk.END)
        for i in range(wf_full.shape[1]):
            self.scintillation_listbox.insert(tk.END, f"Channel {i+1}")
        self.scintillation_listbox.select_set(0)

    # Clear memory (only will run when event/run is switched)
    def reset_analysis_cache(self):
        self.analysis_cache = None
        self.analyzed_triggers = None
        self.loaded_trigger_range = None
        self.trigger_range_start_var.set(0)
        # Hard coded trigger range to be 100 for now
        self.trigger_range_end_var.set(100)

    def new_channel(self):
        # If no channel is selected go to channel 1 and draw fastdaq will clear graphs
        selections = self.scintillation_listbox.curselection()
        if not selections:
            self.draw_fastdaq_scintillation()
            selections = self.scintillation_listbox.curselection() or [0]
            idx = selections[0]
            selections = self.scintillation_listbox.curselection() or [0]
            idx = selections[0]
        idx = selections[0]  # Use first selected index for now
        if idx < 0:
            idx = 0

        # Find how many triggers 
        if 'length' in self.scint_fastdaq_event['scintillation']:
            n_trig = self.scint_fastdaq_event['scintillation']['length']
        else:
            n_trig = self.scint_fastdaq_event['scintillation']['Waveforms'].shape[0]

        n_chan = self.scint_fastdaq_event["scintillation"]["Waveforms"].shape[1]

        # Determine trigger window for analysis
        start_trigger = self.trigger_range_start_var.get()
        end_trigger = self.trigger_range_end_var.get()

        # Clamp to valid range
        start_trigger = max(0, min(start_trigger, n_trig - 1))
        end_trigger = max(start_trigger + 1, min(end_trigger, n_trig))

        self.window_start = start_trigger

        # Check if trigger range exists or new range is different from current
        reload_analysis = (
            not hasattr(self, 'loaded_trigger_range') or 
            self.loaded_trigger_range != (start_trigger, end_trigger)
        )
        if reload_analysis:
            # Make sure waveforms is in the right format (funciton)
            wf_val = self.scint_fastdaq_event["scintillation"]["Waveforms"]
            if callable(wf_val):
                # Set start and end triggers 
                batch_ev = GetScint(
                    self.scint_fastdaq_event,
                    start=start_trigger,
                    end=end_trigger
                )
            else:
                ev = self.scint_fastdaq_event
                batch_ev = { k:v for k,v in ev.items() if k!="scintillation" }
                scint_funcs = {}
                for key, val in ev["scintillation"].items():
                    # Metadata stays as-is
                    if key in ("loaded", "length", "sample_rate", "EventCounter"):
                        scint_funcs[key] = val
                    else:
                        if callable(val):
                            # Already a loader
                            scint_funcs[key] = val
                        else:
                            arr = val
                            # Wrap array in a function
                            scint_funcs[key] = (
                                lambda start, end, length=None, _a=arr: _a[start:end]
                            )
                # Now all scint_funcs[key] are callable
                scint_funcs["length"] = end_trigger - start_trigger
                batch_ev["scintillation"] = scint_funcs
            # Pull cut array
            new_pulses  = SiPMPulsesBatched(batch_ev, nwvf_batch=500)

            if self.analysis_cache is None:
                self.analysis_cache = {
                    key: np.full((n_chan, n_trig), np.nan)
                    for key in new_pulses
                }
                self.analyzed_triggers = np.zeros(n_trig, dtype=bool)
    
            # Merge the entire window into the cache
            idxs = np.arange(start_trigger, end_trigger)
            for key in new_pulses:
                self.analysis_cache[key][:, idxs] = new_pulses[key]
            self.analyzed_triggers[start_trigger:end_trigger] = True
            # Remember what range we’ve just loaded
            self.loaded_trigger_range = (start_trigger, end_trigger)

            self.photon = PhotonT0(self.analysis_cache)

        # Now slice out exactly that window from your cache
        self.pulses = {
            k: self.analysis_cache[k][:, start_trigger:end_trigger]
            for k in self.analysis_cache
        }
    



        # Pull waveforms for a trigger across only selected channels
        waveforms = self.scint_fastdaq_event['scintillation']['Waveforms'][self.trigger_index]
        selected_channels = self.scintillation_listbox.curselection()
        if not selected_channels:
            selected_channels = [0]
        all_selected_data = np.array([waveforms[idx] for idx in selected_channels])
        # Pull sampling rate and triggers
        self.data = self.scint_fastdaq_event['scintillation']['Waveforms'][self.trigger_index][idx]
        self.time = np.arange(len(self.data)) * (1 / self.scint_fastdaq_event['scintillation']['sample_rate'])
        num_trigs = self.scint_fastdaq_event['scintillation']['Waveforms'].shape[0]
        self.trigger_count_label.config(text=f"Triggers: {num_trigs}")
        # Pull voltage and time slider values
        self.dt   = self.time[1] - self.time[0] 
        self.t0   = self.time[0]    
        self.t1 = self.time[-1]
        self.v0 = np.min(all_selected_data)
        self.v1 = np.max(all_selected_data)
        self.dv = (self.v1 - self.v0) / 100.0 
        # Update voltage and time slider range
        self.t_start_slider.config(from_=self.t0, to=self.t1, resolution=self.dt)
        self.t_end_slider.config(from_=self.t0, to=self.t1, resolution=self.dt)
        self.v_lower_slider.config(from_=self.v0, to=self.v1, resolution=self.dv)
        self.v_upper_slider.config(from_=self.v0, to=self.v1, resolution=self.dv)
        # Set sliders to max and min of data set
        self.t_start_slider.set(self.t0)
        self.t_end_slider.set(self.t1)
        # Only update values if lock is not checked
        if not self.lock_voltage_var.get():
            self.v_lower_var.set(self.v0)
            self.v_upper_var.set(self.v1)


        # Update FFT slider range
        nyquist = 1 / (2 * self.dt)
        self.f_low_slider.config(from_=0, to=nyquist, resolution=nyquist / 100)
        self.f_high_slider.config(from_=0, to=nyquist, resolution=nyquist / 100)
        self.f_low_slider.set(0)
        self.f_high_slider.set(nyquist)

        self.draw_fastdaq_scintillation()

    # Plotting
    def draw_fastdaq_scintillation(self, val=None):
        # If no event selected
        if self.scint_fastdaq_event == None:
            self.scint_error()
            return
        # If load SiPM not checked
        if not self.load_fastdaq_scintillation_var.get():
            self.scintillation_tab_right.grid_forget()
            return
        # If event loaded is not the same as selected event
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
        self.scintillation_ax.clear()
        self.fft_ax.clear()
        self.gain_ax.clear()
        selections = self.scintillation_listbox.curselection()
        # If no channels selected clear graphs
        if not selections:
            self.scintillation_ax.clear()
            self.fft_ax.clear()
            self.gain_ax.clear()

            self.scintillation_ax.set_title("No channels selected")
            self.scintillation_ax.set_xlabel('')
            self.scintillation_ax.set_ylabel('')

            self.fft_ax.set_title("No FFT Data")
            self.fft_ax.set_xlabel('')
            self.fft_ax.set_ylabel('')

            self.gain_ax.set_title("No Gain Data")
            self.gain_ax.set_xlabel('')
            self.gain_ax.set_ylabel('')

            self.scintillation_canvas.draw_idle()
            return

        vmins = []
        vmaxs = []
        fft_mins = []
        fft_maxs = []


        for idx in selections:
            data = self.scint_fastdaq_event['scintillation']['Waveforms'][self.trigger_index][idx]
            time = np.arange(len(data)) * (1 / self.scint_fastdaq_event['scintillation']['sample_rate'])

            filtered = self.filter_signal_by_freq(data, self.f_low_var.get(), self.f_high_var.get())
            self.scintillation_ax.plot(time, data, label=f'Raw Ch {idx + 1}')
            self.scintillation_ax.plot(time, filtered, linestyle='--', label=f'Filtered Ch {idx + 1}')

            vmins.append(np.min(data))
            vmaxs.append(np.max(data))

            data = self.scint_fastdaq_event['scintillation']['Waveforms'][self.trigger_index][idx]
            sample_rate = self.scint_fastdaq_event['scintillation']['sample_rate']
            dt = 1 / sample_rate

            filtered = self.filter_signal_by_freq(data, self.f_low_var.get(), self.f_high_var.get())
            fft_vals = np.fft.rfft(filtered)
            freqs = np.fft.rfftfreq(len(filtered), d=dt)
            fft_mag = np.abs(fft_vals)

            self.fft_ax.plot(freqs[1:], fft_mag[1:], label=f'Ch {idx}')
            fft_mins.append(np.min(fft_mag[1:]))
            fft_maxs.append(np.max(fft_mag[1:]))

        #     try:
        #         t0_val = self.pulses['hit_t0'][idx, self.trigger_index]
        #         amp_val = self.pulses['hit_amp'][idx, self.trigger_index]
        #         if not np.isnan(t0_val) and amp_val > 0:
        #             self.scintillation_ax.axvline(t0_val, color='red', linestyle='dotted', label=f'hit_t0 Ch {idx}')
        #     except Exception as e:
        #         print(f"Failed to add hit_t0 for channel {idx}: {e}")
            # try:
            #     baseline_val = self.pulses['baseline'][idx, self.trigger_index]
            #     rms_val = self.pulses['rms'][idx, self.trigger_index]
            #     threshold_val = baseline_val + 5 * rms_val
            #     self.scintillation_ax.axhline(threshold_val, color='purple', linestyle='dotted', label=f'Thresh Ch {idx}')
            # except Exception as e:
            #     print(f"Failed to add threshold line for Ch {idx}: {e}")
        # Waveform axis settings
        self.scintillation_ax.relim()
        self.scintillation_ax.autoscale_view()
        self.scintillation_ax.set_xlim(start, end)
        self.scintillation_ax.set_ylim(vlow, vhigh)
        self.scintillation_ax.set_title("Channels: " + ", ".join(str(i+1) for i in selections) + " Run: " + str(self.run) + " Event: " + str(self.event))
        self.scintillation_ax.set_xlabel('[s]')
        self.scintillation_ax.set_ylabel('[mV]')
        self.scintillation_ax.legend(loc='upper right', fontsize='small')

        # Histogram of hit amplitudes
        amps = self.photon['amp']
        mask = getattr(self, 'analyzed_triggers', None)
        if mask is None:
            # fallback to nan‐filter only
            valid_amps = amps[~np.isnan(amps)]
        else:
            # only include the triggers we've loaded into the cache
            loaded_amps = amps[mask]
            valid_amps  = loaded_amps[~np.isnan(loaded_amps)]

        self.gain_ax.hist(valid_amps, bins=50)
        self.gain_ax.set_title(f"Hits per Amplitude histogram (n={len(valid_amps)})")
        mask = getattr(self, 'analyzed_triggers', None)
        if mask is None:
            # fallback to nan‐filter only
            valid_amps = amps[~np.isnan(amps)]
        else:
            # only include the triggers we've loaded into the cache
            loaded_amps = amps[mask]
            valid_amps  = loaded_amps[~np.isnan(loaded_amps)]

        self.gain_ax.hist(valid_amps, bins=50)
        self.gain_ax.set_title(f"Hits per Amplitude histogram (n={len(valid_amps)})")
        self.gain_ax.set_xlabel("Pulse amplitude (mV)")
        self.gain_ax.set_ylabel("Hits")

        # Plot FFT
        if fft_mins and fft_maxs:
            self.fft_ax.set_ylim(min(fft_mins), max(fft_maxs))

        self.fft_ax.set_xlim(flow, fhigh)
        self.fft_ax.set_title("FFT Magnitude")
        self.fft_ax.set_xlabel("Frequency (Hz)")
        self.fft_ax.set_ylabel("Magnitude")

        if selections:
            self.fft_ax.legend(loc='upper right', fontsize='small')

                # Show updated channel variable info and cuts
        self.update_channel_info_display()

        # Render
        self.scintillation_canvas.draw_idle()
        self.scintillation_canvas.get_tk_widget().grid(row=0, column=1, sticky='NW')

    # Error handling
    def scint_error(self):
        # Clear event and combobox
        self.scint_fastdaq_event = None
        self.scintillation_listbox.delete(0, tk.END)
        # Clear graphs and display error message with run and event for bug testing
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

    # Get FFT arrays
    def filter_signal_by_freq(self, data, flow, fhigh):
        fft = np.fft.rfft(data)
        freqs = np.fft.rfftfreq(len(data), d=self.dt)

        # Zero out frequencies outside the desired range
        bandpass = (freqs >= flow) & (freqs <= fhigh)
        fft[~bandpass] = 0

        # Inverse FFT to return filtered time signal
        return np.fft.irfft(fft, n=len(data))

    # Wrapper function for sliders
    def on_slider_release(self, var):
        self.draw_fastdaq_scintillation()   

    # If trig entry is below or above the range of trig
    def on_trigger_entry_change(self, event):
        try:
            idx = int(self.trigger_var.get()) - 1
            max_idx = self.scint_fastdaq_event['scintillation']['Waveforms'].shape[0]
            if 0 <= idx < max_idx:
                self.trigger_index = idx
                self.new_channel()
            else:
                print(f"Trigger index {idx + 1} out of range.")
        except ValueError:
            print("Invalid trigger index entered.")

    # Shifting triggers by button value
    def shift_trigger(self, step):
        max_idx = self.scint_fastdaq_event['scintillation']['Waveforms'].shape[0]
        new_idx = self.trigger_index + step
        new_idx = max(0, min(new_idx, max_idx - 1))
        self.trigger_index = new_idx
        self.trigger_var.set(str(new_idx + 1))
        self.new_channel()
    
    def update_channel_info_display(self):
        for widget in self.info_frame.winfo_children():
            widget.destroy()

        selected_channels = self.scintillation_listbox.curselection()
        if not selected_channels:
            return

        # Variable names on left
        variable_names = [
            "hit_t0", "hit_amp", "hit_area", "wvf_area",
            "Second Pulse", "Baseline", "RMS"
        ]
        for row, var in enumerate(variable_names, start=1):
            tk.Label(self.info_frame, text=var).grid(row=row, column=0, sticky='w')

        # Fill in channel headers + values
        for col, ch_idx in enumerate(selected_channels, start=1):
            tk.Label(self.info_frame, text=f"Ch {ch_idx + 1}", font=("Arial", 9, "underline")).grid(row=0, column=col, padx=5)
            try:
                values = [
                    f"{self.pulses['hit_t0'][ch_idx, self.trigger_index]:.4e}",
                    f"{self.pulses['hit_amp'][ch_idx, self.trigger_index]:.4e}",
                    f"{self.pulses['hit_area'][ch_idx, self.trigger_index]:.4e}",
                    f"{self.pulses['wvf_area'][ch_idx, self.trigger_index]:.4e}",
                    "Yes" if self.pulses['second_pulse'][ch_idx, self.trigger_index] else "No",
                    f"{self.pulses['baseline'][ch_idx, self.trigger_index]:.4f}",
                    f"{self.pulses['rms'][ch_idx, self.trigger_index]:.4f}",
                ]
                for row, val in enumerate(values, start=1):
                    tk.Label(self.info_frame, text=val).grid(row=row, column=col, padx=5)
            except Exception as e:
                print(f"Error showing pulse info for Channel {ch_idx}: {e}")


        self.cut_start_row = row  # store where to begin the cut widgets
        # Recreate cut widgets below channel info
        self.cut_widgets = {}

        cut_variables = ["hit_amp", "hit_area", "wvf_area", "hit_t0", "second_pulse"]
        for i, var in enumerate(cut_variables):
            row = self.cut_start_row + i
            if var == "second_pulse":
                yes_var = tk.BooleanVar(value=True)
                no_var = tk.BooleanVar(value=True)
                yes_check = ttk.Checkbutton(self.info_frame, text="Yes", variable=yes_var)
                no_check = ttk.Checkbutton(self.info_frame, text="No", variable=no_var)

                if self.enable_cuts_var.get():
                    yes_check.grid(row=row, column=3)
                    no_check.grid(row=row, column=4)

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

                if self.enable_cuts_var.get():
                    entry.grid(row=row, column=3)
                    above_button.grid(row=row, column=4)
                    below_button.grid(row=row, column=5)

                self.cut_widgets[var] = {
                    "entry": entry,
                    "mode": mode,
                    "above": above_button,
                    "below": below_button
                }



    # Creating figure, canvas, and axes for plots
    def scintillation_canvas_setup(self):
        # Figure and canvas for plotting
        self.scintillation_fig = Figure(figsize=(7, 9), dpi=100)
        gs = self.scintillation_fig.add_gridspec(3, 1)

        # Creating axes in fig frame
        self.scintillation_ax = self.scintillation_fig.add_subplot(gs[0, 0])
        self.gain_ax = self.scintillation_fig.add_subplot(gs[2, 0])
        self.fft_ax = self.scintillation_fig.add_subplot(gs[1, 0])

        # Spacing and fitting
        self.scintillation_fig.subplots_adjust(left=0.12, right=0.98, top=0.95, bottom=0.05, hspace=0.8)

        # Creating a frame with plots. If you want multiple plots side by side create a new figure and put axes on it and call as so
        self.scintillation_canvas = FigureCanvasTkAgg(self.scintillation_fig, self.scintillation_tab_right)
        self.scintillation_canvas.get_tk_widget().grid(row=0, column=1, sticky='NSEW')

        self.scintillation_tab_right.grid_rowconfigure(0, weight=1)
        self.scintillation_tab_right.grid_columnconfigure(1, weight=1)

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
        self.info_frame.grid(row=7, column=0, columnspan=4, sticky='WE', pady=(10, 0))
        self.info_frame.grid(row=7, column=0, columnspan=4, sticky='WE', pady=(10, 0))
        self.trigger_step_frame = tk.Frame(self.scintillation_tab_left)
        self.trigger_step_frame.grid(row=1, column=4, padx=(10, 0), sticky="W")

        # Load SiPM checkbutton
        self.load_sipm_checkbutton = tk.Checkbutton(
            self.scintillation_tab_left,

            text='Load SiPM',

            variable=self.load_fastdaq_scintillation_var,
            command=self.load_fastdaq_scintillation
        )
        self.load_sipm_checkbutton.grid(row=0, column=0, columnspan=2, sticky='WE')

        # Channel selector combobox
        tk.Label(self.scintillation_tab_left, text='Channel:').grid(row=1, column=0, sticky='WE')
        # Frame to hold Listbox + Scrollbar
        listbox_frame = tk.Frame(self.scintillation_tab_left)
        listbox_frame.grid(row=1, column=1, sticky='WE')
        # Listbox
        self.scintillation_listbox = tk.Listbox(
            listbox_frame,
            selectmode='multiple',
            height=6,
            exportselection=False,
            yscrollcommand=lambda *args: self.listbox_scrollbar.set(*args)
        )
        self.scintillation_listbox.pack(side=tk.LEFT, fill=tk.BOTH)
        # Scrollbar
        self.listbox_scrollbar = tk.Scrollbar(
            listbox_frame,
            orient="vertical",
            command=self.scintillation_listbox.yview
        )
        self.listbox_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.scintillation_listbox.bind("<<ListboxSelect>>", lambda _: self.new_channel())
    
        # Label showing number of triggers
        self.trigger_count_label = tk.Label(self.scintillation_tab_left, text="Triggers: ?")
        self.trigger_count_label.grid(row=1, column=2, padx=(10, 0), sticky="W")

        # Entry box for trigger selection
        self.trigger_var = tk.StringVar()
        self.trigger_var.set("1")
        self.trigger_entry = tk.Entry(self.scintillation_tab_left, textvariable=self.trigger_var, width=6)
        self.trigger_entry.grid(row=1, column=3, padx=(2, 10), sticky="W")
        self.trigger_entry.bind("<Return>", self.on_trigger_entry_change)

        # 6 skip buttons in 3x2 grid for triggers
        btn_specs = [
            ("-1", -1),   ("+1", 1),
            ("-10", -10), ("+10", 10),
            ("-100", -100), ("+100", 100),
            ("-1000", -1000), ("+1000", 1000),
            ("-10000", -10000), ("+10000", -10000)
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

        # Cut widgets
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

        # Entry fields for trigger range (start and end)
        tk.Label(self.scintillation_tab_left, text='Start Trigger:').grid(row=2, column=0, sticky='E')
        self.trigger_range_start_var = tk.IntVar(value=0)
        self.trigger_range_start_entry = tk.Entry(self.scintillation_tab_left, textvariable=self.trigger_range_start_var, width=6)
        self.trigger_range_start_entry.grid(row=2, column=1, sticky='W')

        tk.Label(self.scintillation_tab_left, text='End Trigger:').grid(row=2, column=2, sticky='E')
        self.trigger_range_end_var = tk.IntVar(value=0)
        self.trigger_range_end_entry = tk.Entry(self.scintillation_tab_left, textvariable=self.trigger_range_end_var, width=6)
        self.trigger_range_end_entry.grid(row=2, column=3, sticky='W')

        # Entry fields for trigger range (start and end)
        tk.Label(self.scintillation_tab_left, text='Start Trigger:').grid(row=2, column=0, sticky='E')
        self.trigger_range_start_var = tk.IntVar(value=0)
        self.trigger_range_start_entry = tk.Entry(self.scintillation_tab_left, textvariable=self.trigger_range_start_var, width=6)
        self.trigger_range_start_entry.grid(row=2, column=1, sticky='W')

        tk.Label(self.scintillation_tab_left, text='End Trigger:').grid(row=2, column=2, sticky='E')
        self.trigger_range_end_var = tk.IntVar(value=0)
        self.trigger_range_end_entry = tk.Entry(self.scintillation_tab_left, textvariable=self.trigger_range_end_var, width=6)
        self.trigger_range_end_entry.grid(row=2, column=3, sticky='W')

        # Time domain slider
        self.t_start_var = tk.DoubleVar(value=0.0)
        self.t_end_var   = tk.DoubleVar(value=0.0)

        tk.Label(self.scintillation_tab_left, text='T start:').grid(row=3, column=0, sticky='E')
        tk.Label(self.scintillation_tab_left, text='T start:').grid(row=3, column=0, sticky='E')
        self.t_start_slider = tk.Scale(
            self.scintillation_tab_left,
            variable=self.t_start_var,
            orient='horizontal',
            from_=0, to=0,         
            resolution=1,             
            showvalue=True,
            length=200
        )
        self.t_start_slider.grid(row=3, column=1, sticky='WE')
        self.t_start_slider.grid(row=3, column=1, sticky='WE')
        self.t_start_slider.bind("<ButtonRelease-1>", self.on_slider_release)


        tk.Label(self.scintillation_tab_left, text='T end:').grid(row=3, column=2, sticky='E')
        tk.Label(self.scintillation_tab_left, text='T end:').grid(row=3, column=2, sticky='E')
        self.t_end_slider = tk.Scale(
            self.scintillation_tab_left,
            variable=self.t_end_var,
            orient='horizontal',
            from_=0, to=0,            
            resolution=1,            
            showvalue=True,
            length=200
        )
        self.t_end_slider.grid(row=3, column=3, sticky='WE')
        self.t_end_slider.grid(row=3, column=3, sticky='WE')
        self.t_end_slider.bind("<ButtonRelease-1>", self.on_slider_release)


        # Voltage range sliders
        self.v_lower_var = tk.DoubleVar(value=0.0)
        self.v_upper_var = tk.DoubleVar(value=0.0)

        tk.Label(self.scintillation_tab_left, text='V lower:').grid(row=4, column=0, sticky='E')
        tk.Label(self.scintillation_tab_left, text='V lower:').grid(row=4, column=0, sticky='E')
        self.v_lower_slider = tk.Scale(
            self.scintillation_tab_left,
            variable=self.v_lower_var,
            orient='horizontal',
            from_=0, to=0,         
            resolution=1e-3,        
            showvalue=True,
            length=200
        )
        self.v_lower_slider.grid(row=4, column=1, sticky='WE')
        self.v_lower_slider.grid(row=4, column=1, sticky='WE')
        self.v_lower_slider.bind("<ButtonRelease-1>", self.on_slider_release)

        tk.Label(self.scintillation_tab_left, text='V upper:').grid(row=4, column=2, sticky='E')
        tk.Label(self.scintillation_tab_left, text='V upper:').grid(row=4, column=2, sticky='E')
        self.v_upper_slider = tk.Scale(
            self.scintillation_tab_left,
            variable=self.v_upper_var,
            orient='horizontal',
            from_=0, to=0,         
            resolution=1e-3,    
            showvalue=True,
            length=200
        )
        self.v_upper_slider.grid(row=4, column=3, sticky='WE')
        self.v_upper_slider.grid(row=4, column=3, sticky='WE')
        self.v_upper_slider.bind("<ButtonRelease-1>", self.on_slider_release)

        self.lock_voltage_var = tk.BooleanVar(value=False)
        self.lock_voltage_check = tk.Checkbutton(
            self.scintillation_tab_left,
            text='Lock Voltage',
            variable=self.lock_voltage_var
        )
        self.lock_voltage_check.grid(row=4, column=4, padx=(10, 0), sticky='W')
        self.lock_voltage_check.grid(row=4, column=4, padx=(10, 0), sticky='W')


        # Frequency cutoff sliders
        self.f_low_var = tk.DoubleVar(value=0.0)
        self.f_high_var = tk.DoubleVar(value=0.0)

        tk.Label(self.scintillation_tab_left, text='F low (Hz):').grid(row=5, column=0, sticky='E')
        tk.Label(self.scintillation_tab_left, text='F low (Hz):').grid(row=5, column=0, sticky='E')
        self.f_low_slider = tk.Scale(
            self.scintillation_tab_left,
            variable=self.f_low_var,
            orient='horizontal',
            from_=0, to=0,
            resolution=1,
            showvalue=True,
            length=200
        )
        self.f_low_slider.grid(row=5, column=1, sticky='WE')
        self.f_low_slider.grid(row=5, column=1, sticky='WE')
        self.f_low_slider.bind("<ButtonRelease-1>", self.on_slider_release)

        tk.Label(self.scintillation_tab_left, text='F high (Hz):').grid(row=5, column=2, sticky='E')
        tk.Label(self.scintillation_tab_left, text='F high (Hz):').grid(row=5, column=2, sticky='E')
        self.f_high_slider = tk.Scale(
            self.scintillation_tab_left,
            variable=self.f_high_var,
            orient='horizontal',
            from_=0, to=0,
            resolution=1,
            showvalue=True,
            length=200
        )
        self.f_high_slider.grid(row=5, column=3, sticky='WE')
        self.f_high_slider.grid(row=5, column=3, sticky='WE')
        self.f_high_slider.bind("<ButtonRelease-1>", self.on_slider_release)

        # Reload button
        self.reload_fastdaq_scintillation_button = tk.Button(
            self.scintillation_tab_left,
            text='Reload',
            command=self.new_channel
        )
        self.reload_fastdaq_scintillation_button.grid(row=6, column=0, columnspan=4, sticky='WE')
        self.reload_fastdaq_scintillation_button.grid(row=6, column=0, columnspan=4, sticky='WE')

    # Cut widget show/hide
    def place_cut_widgets(self):
        if self.enable_cuts_var.get():
            for i, (var, widgets) in enumerate(self.cut_widgets.items()):
                row = self.cut_start_row + i
                if var == "second_pulse":
                    widgets["yes_check"].grid(row=row, column=3)
                    widgets["no_check"].grid(row=row, column=4)
                else:
                    widgets["entry"].grid(row=row, column=3)
                    widgets["above"].grid(row=row, column=4)
                    widgets["below"].grid(row=row, column=5)
        else:
            for widgets in self.cut_widgets.values():
                if "yes_check" in widgets:
                    widgets["yes_check"].grid_remove()
                    widgets["no_check"].grid_remove()
                else:
                    widgets["entry"].grid_remove()
                    widgets["above"].grid_remove()
                    widgets["below"].grid_remove()

    def toggle_cut_widgets(self):
        self.update_channel_info_display()
