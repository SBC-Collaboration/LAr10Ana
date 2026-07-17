import gc
import os
import sys
import numpy as np
import tkinter as tk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Hacky
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
from GetEvent import GetEvent, GetScint
from sbcbinaryformat import Streamer

# Recon scintillation.sbc columns surfaced in the Pulse Info panel. Each is stored
# per channel (shape (32,)) per CAEN trigger row.
INFO_FIELDS = ["hit_t0", "hit_amp", "hit_area", "wvf_area", "second_pulse", "baseline", "rms"]


class Scintillation(tk.Frame):
    """SiPM tab.

    Waveform/FFT (top two axes) are lazy read on demand from the *raw*
    scintillation.sbc(not recon file) one trigger at a time
    The Pulse Info panel and amplitude histogram are
    read from the precomputed recon scintillation.sbc
    """

    def __init__(self, master=None):
        # State variable for loading SiPM data
        self.load_fastdaq_scintillation_var = tk.BooleanVar(value=False)

        # Displayed trigger (0-based)
        self.trigger_index = 0

        # Raw event (lazy) used only to draw the selected trigger's waveform/FFT
        self.scint_fastdaq_event = None
        self.sample_rate = None
        self.dt = 1.0
        self.n_trigs = 0
        self.current_wf = None  # (n_chs, rec_len) for the displayed trigger

        # Precomputed recon analysis is one binary per run.
        # self._recon_all is the full run, self.recon is the current event.
        self._recon_all = None
        self.recon_run = None
        self.recon = None

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

        self.trigger_index = 0
        self.trigger_var.set("1")
        try:
            self.load_recon()
            self.load_event()
            if self.scint_fastdaq_event is None and self.recon is None:
                self.scint_error()
                return
            self.populate_channel_listbox()
            self.refresh()
        except Exception as e:
            print(f"Scintillation load failed for run {self.run} event {self.event}: {e}")
            self.scint_error()
        gc.collect()

    def load_recon(self):
        # Load the precomputed recon scintillation.sbc for this run (cached), then filter
        # to the current event. Rows = CAEN triggers; columns 0-31 = channels.
        if self.recon_run != self.run or self._recon_all is None:
            self._recon_all = None
            self.recon_run = self.run
            path = self._find_recon('scintillation.sbc', self.run)
            if path is not None:
                self._recon_all = Streamer(path).data

        if self._recon_all is None:
            self.recon = None
            return

        d = self._recon_all
        if 'ev' in d.dtype.names:
            self.recon = d[d['ev'] == int(self.event)]
        else:
            self.recon = d

    def load_event(self):
        # Lazy-load the raw event so large runs don't pull every waveform into memory.
        # run_control is needed for sample_rate; event_info kept for parity with GetEvent.
        run_path = os.path.join(self.raw_directory, self.run)
        selected = ["run_control", "scintillation", "event_info"]
        event = GetEvent(run_path, self.event, *selected,
                         strictMode=False, lazy_load_scintillation=True)
        scint = event["scintillation"]
        if not scint.get("loaded"):
            self.scint_fastdaq_event = None
            self.sample_rate = None
            self.n_trigs = 0
            return
        self.scint_fastdaq_event = event
        self.sample_rate = scint.get("sample_rate")
        self.n_trigs = int(scint.get("length", 0))
        # Seed the histogram trigger window to cover the whole event
        self.trigger_range_start_var.set(0)
        self.trigger_range_end_var.set(self.n_trigs)

    def get_waveforms(self, trig):
        # Read a single trigger's raw waveforms: (n_chs, rec_len), or None.
        if self.scint_fastdaq_event is None or self.n_trigs == 0:
            return None
        trig = max(0, min(trig, self.n_trigs - 1))
        one = GetScint(self.scint_fastdaq_event, start=trig, end=trig + 1)
        wf = np.asarray(one['scintillation']['Waveforms'])
        return wf[0] if wf.ndim == 3 else wf

    def populate_channel_listbox(self):
        # Channel count from the raw waveform, falling back to the recon column count.
        wf = self.get_waveforms(self.trigger_index)
        if wf is not None:
            n_chan = wf.shape[0]
        elif self.recon is not None and len(self.recon):
            n_chan = np.asarray(self.recon['hit_amp']).shape[1]
        else:
            n_chan = 0
        self.scintillation_listbox.delete(0, tk.END)
        for i in range(n_chan):
            self.scintillation_listbox.insert(tk.END, f"Channel {i+1}")
        if n_chan:
            self.scintillation_listbox.select_set(0)

    def refresh(self, *_):
        # Fetch the current trigger's waveform, update slider ranges, and redraw.
        self.current_wf = self.get_waveforms(self.trigger_index)
        if self.sample_rate:
            self.dt = 1.0 / self.sample_rate
        self.trigger_count_label.config(text=f"Triggers: {self.n_trigs}")
        self.update_waveform_settings()
        self.draw_fastdaq_scintillation()

    def update_waveform_settings(self):
        # Configure the time/voltage/FFT slider ranges from the displayed trigger.
        wf = self.current_wf
        if wf is None:
            return
        selected = self.scintillation_listbox.curselection() or [0]
        selected = [c for c in selected if c < wf.shape[0]] or [0]
        sel_data = wf[selected]

        rec_len = wf.shape[1]
        time = np.arange(rec_len) * self.dt
        self.t0, self.t1 = time[0], time[-1]
        dt = self.dt if self.dt else 1.0

        self.v0 = float(np.min(sel_data))
        self.v1 = float(np.max(sel_data))
        dv = (self.v1 - self.v0) / 100.0 or 1e-3

        self.t_start_slider.config(from_=self.t0, to=self.t1, resolution=dt)
        self.t_end_slider.config(from_=self.t0, to=self.t1, resolution=dt)
        self.v_lower_slider.config(from_=self.v0, to=self.v1, resolution=dv)
        self.v_upper_slider.config(from_=self.v0, to=self.v1, resolution=dv)
        self.t_start_slider.set(self.t0)
        self.t_end_slider.set(self.t1)
        if not self.lock_voltage_var.get():
            self.v_lower_var.set(self.v0)
            self.v_upper_var.set(self.v1)

        nyquist = 1.0 / (2 * dt)
        self.f_low_slider.config(from_=0, to=nyquist, resolution=nyquist / 100)
        self.f_high_slider.config(from_=0, to=nyquist, resolution=nyquist / 100)
        self.f_low_slider.set(0)
        self.f_high_slider.set(nyquist)

    def draw_fastdaq_scintillation(self, *_):
        if not self.load_fastdaq_scintillation_var.get():
            self.scintillation_tab_right.grid_forget()
            return
        self.scintillation_tab_right.grid(row=0, column=1, sticky='NW')

        self.scintillation_ax.clear()
        self.fft_ax.clear()
        self.gain_ax.clear()

        selections = self.scintillation_listbox.curselection()
        wf = self.current_wf

        # Waveform + FFT(raw scintillation.sbc for the selected trigger and channels)
        if wf is not None and selections:
            flow = self.f_low_var.get()
            fhigh = self.f_high_var.get()
            fft_mins, fft_maxs = [], []
            for idx in selections:
                if idx >= wf.shape[0]:
                    continue
                data = wf[idx].astype(float)
                time = np.arange(len(data)) * self.dt
                filtered = self.filter_signal_by_freq(data, flow, fhigh)
                self.scintillation_ax.plot(time, data, label=f'Raw Ch {idx + 1}')
                self.scintillation_ax.plot(time, filtered, linestyle='--', label=f'Filtered Ch {idx + 1}')

                fft_mag = np.abs(np.fft.rfft(filtered))
                freqs = np.fft.rfftfreq(len(filtered), d=self.dt)
                self.fft_ax.plot(freqs[1:], fft_mag[1:], label=f'Ch {idx + 1}')
                fft_mins.append(float(np.min(fft_mag[1:])))
                fft_maxs.append(float(np.max(fft_mag[1:])))

            self.scintillation_ax.set_xlim(self.t_start_var.get(), self.t_end_var.get())
            self.scintillation_ax.set_ylim(self.v_lower_var.get(), self.v_upper_var.get())
            self.scintillation_ax.set_title(
                "Channels: " + ", ".join(str(i + 1) for i in selections)
                + f"  Run: {self.run}  Event: {self.event}  Trigger: {self.trigger_index + 1}")
            self.scintillation_ax.set_xlabel('[s]')
            self.scintillation_ax.set_ylabel('[ADC]')
            self.scintillation_ax.legend(loc='upper right', fontsize='small')

            if fft_mins and fft_maxs:
                self.fft_ax.set_ylim(min(fft_mins), max(fft_maxs))
            self.fft_ax.set_xlim(flow, fhigh)
            self.fft_ax.set_title("FFT Magnitude")
            self.fft_ax.set_xlabel("Frequency (Hz)")
            self.fft_ax.set_ylabel("Magnitude")
            self.fft_ax.legend(loc='upper right', fontsize='small')
        else:
            msg = "No channels selected" if wf is not None else "No raw waveform data"
            self.scintillation_ax.set_title(msg)
            self.fft_ax.set_title("No FFT Data")

        # Amplitude histogram(recon file)
        amps = self.event_amps()
        if amps is None:
            self.gain_ax.set_title("No recon analysis for this run")
        else:
            self.gain_ax.hist(amps, bins=100)
            self.gain_ax.set_title(f"Hits per Amplitude histogram (n={len(amps)})")
            self.gain_ax.set_xlabel("Pulse amplitude (mV)")
            self.gain_ax.set_ylabel("Hits")

        self.update_channel_info_display()
        self.scintillation_canvas.draw_idle()

    def event_amps(self):
        # Flatten recon hit_amp over the selected channels and the trigger window.
        if self.recon is None or len(self.recon) == 0:
            return None
        amps = np.asarray(self.recon['hit_amp'])  # (n_trig_event, n_chan)
        start = max(0, self.trigger_range_start_var.get())
        end = self.trigger_range_end_var.get()
        if end <= start:
            end = amps.shape[0]
        amps = amps[start:min(end, amps.shape[0])]
        sel = self.scintillation_listbox.curselection()
        cols = [c for c in sel if c < amps.shape[1]] or list(range(amps.shape[1]))
        vals = amps[:, cols].ravel()
        return vals[~np.isnan(vals)]

    def scint_error(self):
        self.scint_fastdaq_event = None
        self.recon = None
        self.current_wf = None
        self.scintillation_listbox.delete(0, tk.END)
        for ax, msg in ((self.scintillation_ax, f"No Data For Run {self.run} Event {self.event}"),
                        (self.gain_ax, "No Gain Data"),
                        (self.fft_ax, "No FFT Data")):
            ax.clear()
            ax.text(0.5, 0.5, msg, transform=ax.transAxes, fontsize=10, ha='center', va='center')
        self.scintillation_canvas.draw_idle()

    # Get FFT arrays
    def filter_signal_by_freq(self, data, flow, fhigh):
        fft = np.fft.rfft(data)
        freqs = np.fft.rfftfreq(len(data), d=self.dt)
        bandpass = (freqs >= flow) & (freqs <= fhigh)
        fft[~bandpass] = 0
        return np.fft.irfft(fft, n=len(data))

    # Wrapper function for sliders
    def on_slider_release(self, var):
        self.draw_fastdaq_scintillation()

    def on_trigger_entry_change(self, event):
        try:
            idx = int(self.trigger_var.get()) - 1
        except ValueError:
            print("Invalid trigger index entered.")
            return
        if 0 <= idx < self.n_trigs:
            self.trigger_index = idx
            self.refresh()
        else:
            print(f"Trigger index {idx + 1} out of range.")

    def shift_trigger(self, step):
        if self.n_trigs == 0:
            return
        new_idx = max(0, min(self.trigger_index + step, self.n_trigs - 1))
        self.trigger_index = new_idx
        self.trigger_var.set(str(new_idx + 1))
        self.refresh()

    def update_channel_info_display(self):
        for widget in self.info_frame.winfo_children():
            widget.destroy()

        selected_channels = self.scintillation_listbox.curselection()
        if not selected_channels:
            return

        variable_names = ["hit_t0", "hit_amp", "hit_area", "wvf_area",
                          "Second Pulse", "Baseline", "RMS"]
        for row, var in enumerate(variable_names, start=1):
            tk.Label(self.info_frame, text=var).grid(row=row, column=0, sticky='w')

        have_recon = self.recon is not None and self.trigger_index < len(self.recon)
        for col, ch_idx in enumerate(selected_channels, start=1):
            tk.Label(self.info_frame, text=f"Ch {ch_idx + 1}",
                     font=("Arial", 9, "underline")).grid(row=0, column=col, padx=5)
            if not have_recon:
                tk.Label(self.info_frame, text="—").grid(row=1, column=col, padx=5)
                continue
            r = self.recon[self.trigger_index]
            try:
                values = [
                    f"{r['hit_t0'][ch_idx]:.4e}",
                    f"{r['hit_amp'][ch_idx]:.4e}",
                    f"{r['hit_area'][ch_idx]:.4e}",
                    f"{r['wvf_area'][ch_idx]:.4e}",
                    "Yes" if r['second_pulse'][ch_idx] else "No",
                    f"{r['baseline'][ch_idx]:.4f}",
                    f"{r['rms'][ch_idx]:.4f}",
                ]
                for row, val in enumerate(values, start=1):
                    tk.Label(self.info_frame, text=val).grid(row=row, column=col, padx=5)
            except Exception as e:
                print(f"Error showing pulse info for Channel {ch_idx}: {e}")

    def scintillation_canvas_setup(self):
        self.scintillation_fig = Figure(figsize=(7, 9), dpi=100)
        gs = self.scintillation_fig.add_gridspec(3, 1)
        self.scintillation_ax = self.scintillation_fig.add_subplot(gs[0, 0])
        self.fft_ax = self.scintillation_fig.add_subplot(gs[1, 0])
        self.gain_ax = self.scintillation_fig.add_subplot(gs[2, 0])
        self.scintillation_fig.subplots_adjust(left=0.12, right=0.98, top=0.95, bottom=0.05, hspace=0.8)

        self.scintillation_canvas = FigureCanvasTkAgg(self.scintillation_fig, self.scintillation_tab_right)
        self.scintillation_canvas.get_tk_widget().grid(row=0, column=1, sticky='NSEW')
        self.scintillation_tab_right.grid_rowconfigure(0, weight=1)
        self.scintillation_tab_right.grid_columnconfigure(1, weight=1)

    def create_scintillation_widgets(self):
        # Main frame
        self.scintillation_tab = tk.Frame(self.notebook)
        self.notebook.add(self.scintillation_tab, text='SiPM')

        # Left (controls) and right (plots) panels
        self.scintillation_tab_left = tk.Frame(self.scintillation_tab, bd=5, relief=tk.SUNKEN)
        self.scintillation_tab_left.grid(row=0, column=0, sticky='NW')
        self.scintillation_tab_right = tk.Frame(self.scintillation_tab, bd=5, relief=tk.SUNKEN)
        self.scintillation_tab_right.grid(row=0, column=1, sticky='NW')
        self.info_frame = tk.LabelFrame(self.scintillation_tab_left, text="Pulse Info", padx=5, pady=5)
        self.info_frame.grid(row=7, column=0, columnspan=4, sticky='WE', pady=(10, 0))
        self.trigger_step_frame = tk.Frame(self.scintillation_tab_left)
        self.trigger_step_frame.grid(row=1, column=4, padx=(10, 0), sticky="W")

        # Load SiPM checkbutton
        self.load_sipm_checkbutton = tk.Checkbutton(
            self.scintillation_tab_left, text='Load SiPM',
            variable=self.load_fastdaq_scintillation_var,
            command=self.load_fastdaq_scintillation)
        self.load_sipm_checkbutton.grid(row=0, column=0, columnspan=2, sticky='WE')

        # Channel selector listbox + scrollbar
        tk.Label(self.scintillation_tab_left, text='Channel:').grid(row=1, column=0, sticky='WE')
        listbox_frame = tk.Frame(self.scintillation_tab_left)
        listbox_frame.grid(row=1, column=1, sticky='WE')
        self.scintillation_listbox = tk.Listbox(
            listbox_frame, selectmode='multiple', height=6, exportselection=False,
            yscrollcommand=lambda *args: self.listbox_scrollbar.set(*args))
        self.scintillation_listbox.pack(side=tk.LEFT, fill=tk.BOTH)
        self.listbox_scrollbar = tk.Scrollbar(
            listbox_frame, orient="vertical", command=self.scintillation_listbox.yview)
        self.listbox_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.scintillation_listbox.bind("<<ListboxSelect>>", lambda _: self.refresh())

        # Trigger count + selection entry
        self.trigger_count_label = tk.Label(self.scintillation_tab_left, text="Triggers: ?")
        self.trigger_count_label.grid(row=1, column=2, padx=(10, 0), sticky="W")
        self.trigger_var = tk.StringVar(value="1")
        self.trigger_entry = tk.Entry(self.scintillation_tab_left, textvariable=self.trigger_var, width=6)
        self.trigger_entry.grid(row=1, column=3, padx=(2, 10), sticky="W")
        self.trigger_entry.bind("<Return>", self.on_trigger_entry_change)

        # Trigger step buttons (3x2 grid)
        btn_specs = [
            ("-1", -1), ("+1", 1),
            ("-10", -10), ("+10", 10),
            ("-100", -100), ("+100", 100),
            ("-1000", -1000), ("+1000", 1000),
            ("-10000", -10000), ("+10000", 10000),
        ]
        for i, (label, step) in enumerate(btn_specs):
            btn = tk.Button(self.trigger_step_frame, text=label, width=4,
                            command=lambda s=step: self.shift_trigger(s))
            btn.grid(row=i // 2, column=i % 2, padx=1, pady=1)

        # Histogram trigger window (start/end)
        tk.Label(self.scintillation_tab_left, text='Start Trigger:').grid(row=2, column=0, sticky='E')
        self.trigger_range_start_var = tk.IntVar(value=0)
        tk.Entry(self.scintillation_tab_left, textvariable=self.trigger_range_start_var, width=6).grid(row=2, column=1, sticky='W')
        tk.Label(self.scintillation_tab_left, text='End Trigger:').grid(row=2, column=2, sticky='E')
        self.trigger_range_end_var = tk.IntVar(value=0)
        tk.Entry(self.scintillation_tab_left, textvariable=self.trigger_range_end_var, width=6).grid(row=2, column=3, sticky='W')

        # Time-domain sliders
        self.t_start_var = tk.DoubleVar(value=0.0)
        self.t_end_var = tk.DoubleVar(value=0.0)
        tk.Label(self.scintillation_tab_left, text='T start:').grid(row=3, column=0, sticky='E')
        self.t_start_slider = tk.Scale(self.scintillation_tab_left, variable=self.t_start_var,
                                       orient='horizontal', from_=0, to=0, resolution=1,
                                       showvalue=True, length=200)
        self.t_start_slider.grid(row=3, column=1, sticky='WE')
        self.t_start_slider.bind("<ButtonRelease-1>", self.on_slider_release)
        tk.Label(self.scintillation_tab_left, text='T end:').grid(row=3, column=2, sticky='E')
        self.t_end_slider = tk.Scale(self.scintillation_tab_left, variable=self.t_end_var,
                                     orient='horizontal', from_=0, to=0, resolution=1,
                                     showvalue=True, length=200)
        self.t_end_slider.grid(row=3, column=3, sticky='WE')
        self.t_end_slider.bind("<ButtonRelease-1>", self.on_slider_release)

        # Voltage sliders
        self.v_lower_var = tk.DoubleVar(value=0.0)
        self.v_upper_var = tk.DoubleVar(value=0.0)
        tk.Label(self.scintillation_tab_left, text='V lower:').grid(row=4, column=0, sticky='E')
        self.v_lower_slider = tk.Scale(self.scintillation_tab_left, variable=self.v_lower_var,
                                       orient='horizontal', from_=0, to=0, resolution=1e-3,
                                       showvalue=True, length=200)
        self.v_lower_slider.grid(row=4, column=1, sticky='WE')
        self.v_lower_slider.bind("<ButtonRelease-1>", self.on_slider_release)
        tk.Label(self.scintillation_tab_left, text='V upper:').grid(row=4, column=2, sticky='E')
        self.v_upper_slider = tk.Scale(self.scintillation_tab_left, variable=self.v_upper_var,
                                       orient='horizontal', from_=0, to=0, resolution=1e-3,
                                       showvalue=True, length=200)
        self.v_upper_slider.grid(row=4, column=3, sticky='WE')
        self.v_upper_slider.bind("<ButtonRelease-1>", self.on_slider_release)

        self.lock_voltage_var = tk.BooleanVar(value=False)
        self.lock_voltage_check = tk.Checkbutton(self.scintillation_tab_left, text='Lock Voltage',
                                                 variable=self.lock_voltage_var)
        self.lock_voltage_check.grid(row=4, column=4, padx=(10, 0), sticky='W')

        # Frequency-cutoff sliders
        self.f_low_var = tk.DoubleVar(value=0.0)
        self.f_high_var = tk.DoubleVar(value=0.0)
        tk.Label(self.scintillation_tab_left, text='F low (Hz):').grid(row=5, column=0, sticky='E')
        self.f_low_slider = tk.Scale(self.scintillation_tab_left, variable=self.f_low_var,
                                     orient='horizontal', from_=0, to=0, resolution=1,
                                     showvalue=True, length=200)
        self.f_low_slider.grid(row=5, column=1, sticky='WE')
        self.f_low_slider.bind("<ButtonRelease-1>", self.on_slider_release)
        tk.Label(self.scintillation_tab_left, text='F high (Hz):').grid(row=5, column=2, sticky='E')
        self.f_high_slider = tk.Scale(self.scintillation_tab_left, variable=self.f_high_var,
                                      orient='horizontal', from_=0, to=0, resolution=1,
                                      showvalue=True, length=200)
        self.f_high_slider.grid(row=5, column=3, sticky='WE')
        self.f_high_slider.bind("<ButtonRelease-1>", self.on_slider_release)

        # Reload button (re-reads current trigger + re-applies the histogram window)
        self.reload_fastdaq_scintillation_button = tk.Button(
            self.scintillation_tab_left, text='Reload', command=self.refresh)
        self.reload_fastdaq_scintillation_button.grid(row=6, column=0, columnspan=4, sticky='WE')
