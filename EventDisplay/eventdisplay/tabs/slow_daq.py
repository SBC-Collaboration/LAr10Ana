# Imports
import gc
import json
import os
import matplotlib
import scipy.signal
import tkinter as tk
from tkinter import ttk, DISABLED, NORMAL
import numpy as np
import sys
#
matplotlib.use('TkAgg')
# matplotlib.use('Agg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from PIL import Image, ImageTk
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
from GetEvent import GetEvent
from sbcbinaryformat import Streamer

SLOWDAQ_T0_COLORS = {'pressure': 'r', 'bubble': 'g', 'scint': 'b'}

# CAEN scintillation digitizer clock, ticks/s (matches ana/ScintT0.py)
SCINT_SAMPLE_RATE = 125e6


def _unwrap_caen_timestamp(ts, max_ts=2**31):
    # The CAEN trigger clock is 31 bits and rolls over mid-event. Mirrors
    # ana/ScintT0.py so lollipop times match the t0s computed there; kept local
    # to avoid importing that module's sklearn dependency into the GUI.
    ts = np.asarray(ts, dtype=np.int64)
    rollovers = np.diff(ts, axis=-1, prepend=0) < 0
    return ts + np.cumsum(rollovers, axis=-1) * max_ts


class SlowDAQ(tk.Frame):
    def __init__(self):
        self.slowDAQ_event = None

        self.slowDAQ_ymin = None
        self.slowDAQ_ymax = None
        self.slowDAQ_tmin = None
        self.slowDAQ_tmax = None

        self.slowDAQ_t0_vars = {
            key: tk.BooleanVar(value=False) for key in SLOWDAQ_T0_COLORS
        }

        # Scintillation lollipop overlay. Drawn on a twin y axis because pulse
        # area and the slowDAQ sensor share no units.
        self.slowDAQ_scint_var = tk.BooleanVar(value=False)
        self.slowDAQ_scint_ax = None
        self._scint_lolli_run = None   # run whose recon file is cached
        self._scint_lolli_all = None   # recon scintillation.sbc for the whole run
        self._scint_lolli_cache = {}   # ev -> (x_ms, area) or None

        self.create_slowDAQ_widgets()

    def create_slowDAQ_widgets(self):
        # Tab container
        self.slowDAQ_tab = tk.Frame(self.notebook)
        self.notebook.add(self.slowDAQ_tab, text='Slow DAQ')

        self.slowDAQ_tab.grid_columnconfigure(0, weight=0)
        self.slowDAQ_tab.grid_columnconfigure(1, weight=1)
        self.slowDAQ_tab.grid_rowconfigure(0, weight=1)

        # Left: controls
        self.slowDAQ_tab_left = tk.Frame(self.slowDAQ_tab, bd=5, relief=tk.SUNKEN)
        self.slowDAQ_tab_left.grid(row=0, column=0, sticky='NW')
        self.slowDAQ_tab_left.grid_columnconfigure(0, weight=1)
        self.slowDAQ_tab_left.grid_columnconfigure(1, weight=1)

        self.slowDAQ_load_checkbutton_var = tk.BooleanVar(value=False)
        self.slowDAQ_load_checkbutton = tk.Checkbutton(
            self.slowDAQ_tab_left,
            text='Load Slow DAQ',
            variable=self.slowDAQ_load_checkbutton_var,
            command=self.load_slowDAQ
        )
        self.slowDAQ_load_checkbutton.grid(row=0, column=0, columnspan=2, sticky='we')

        tk.Label(self.slowDAQ_tab_left, text='Sensor:').grid(row=1, column=0, sticky='we')
        self.slowDAQ_combobox = ttk.Combobox(self.slowDAQ_tab_left, width=16, state='disabled')
        self.slowDAQ_combobox.grid(row=1, column=1, sticky='we')

        # When user picks a different sensor, redraw
        self.slowDAQ_combobox.bind('<<ComboboxSelected>>', self.draw_slowDAQ)

        # Data cuts
        self.slowDAQ_tmin_var = tk.StringVar(value="")
        self.slowDAQ_tmax_var = tk.StringVar(value="")
        self.slowDAQ_ymin_var   = tk.StringVar(value="")
        self.slowDAQ_ymax_var   = tk.StringVar(value="")

        # Time window
        tk.Label(self.slowDAQ_tab_left, text='t min [ms]:').grid(row=2, column=0, sticky='we')
        self.slowDAQ_tmin_entry = tk.Entry(
            self.slowDAQ_tab_left, width=10, textvariable=self.slowDAQ_tmin_var
        )
        self.slowDAQ_tmin_entry.grid(row=2, column=1, sticky='we')

        tk.Label(self.slowDAQ_tab_left, text='t max [ms]:').grid(row=3, column=0, sticky='we')
        self.slowDAQ_tmax_entry = tk.Entry(
            self.slowDAQ_tab_left, width=10, textvariable=self.slowDAQ_tmax_var
        )
        self.slowDAQ_tmax_entry.grid(row=3, column=1, sticky='we')

        # Amplitude window
        tk.Label(self.slowDAQ_tab_left, text='y min:').grid(row=4, column=0, sticky='we')
        self.slowDAQ_ymin_entry = tk.Entry(
            self.slowDAQ_tab_left, width=10, textvariable=self.slowDAQ_ymin_var
        )
        self.slowDAQ_ymin_entry.grid(row=4, column=1, sticky='we')

        tk.Label(self.slowDAQ_tab_left, text='y max:').grid(row=5, column=0, sticky='we')
        self.slowDAQ_ymax_entry = tk.Entry(
            self.slowDAQ_tab_left, width=10, textvariable=self.slowDAQ_ymax_var
        )
        self.slowDAQ_ymax_entry.grid(row=5, column=1, sticky='we')

        self.slowDAQ_apply_button = tk.Button(
            self.slowDAQ_tab_left,
            text='Apply cuts',
            command=self.apply_slowDAQ_cuts
        )
        self.slowDAQ_apply_button.grid(row=6, column=0, columnspan=2, sticky='we', pady=(5, 0))

        self.slowDAQ_t0_checkbuttons = {}
        row = 7
        for key in SLOWDAQ_T0_COLORS:
            cb = tk.Checkbutton(
                self.slowDAQ_tab_left,
                text=f'Show {key} t0',
                variable=self.slowDAQ_t0_vars[key],
                command=self.draw_slowDAQ,
                state=DISABLED,
            )
            cb.grid(row=row, column=0, columnspan=2, sticky='we')
            self.slowDAQ_t0_checkbuttons[key] = cb
            row += 1

        self.slowDAQ_scint_checkbutton = tk.Checkbutton(
            self.slowDAQ_tab_left,
            text='Show scint pulses',
            variable=self.slowDAQ_scint_var,
            command=self.draw_slowDAQ,
            state=DISABLED,
        )
        self.slowDAQ_scint_checkbutton.grid(row=row, column=0, columnspan=2, sticky='we')
        row += 1

        # Right: plot
        self.slowDAQ_tab_right = tk.Frame(self.slowDAQ_tab, bd=5, relief=tk.SUNKEN)
        self.slowDAQ_tab_right.grid(row=0, column=1, sticky='NW')
        self.slowDAQ_tab_right.grid_rowconfigure(0, weight=1)
        self.slowDAQ_tab_right.grid_columnconfigure(0, weight=1)
        
        self.slowDAQ_fig = Figure(figsize=(7, 5), dpi=100)
        self.slowDAQ_ax = self.slowDAQ_fig.add_subplot(111)

        self.slowDAQ_canvas = FigureCanvasTkAgg(self.slowDAQ_fig, self.slowDAQ_tab_right)
        self.slowDAQ_canvas.get_tk_widget().grid(row=0, column=0, sticky='nsew')

        # Toolbar for navigation
        self.slowDAQ_toolbar = NavigationToolbar2Tk(
            self.slowDAQ_canvas,
            self.slowDAQ_tab_right,
            pack_toolbar=False,
        )
        self.slowDAQ_toolbar.update()
        self.slowDAQ_toolbar.grid(row=1, column=0, sticky='w')




    def load_slowDAQ(self):
        if not self.slowDAQ_load_checkbutton_var.get():
            # clear plot and combobox
            self.slowDAQ_combobox['values'] = []
            self.slowDAQ_combobox.set('')
            self.slowDAQ_combobox.state(['disabled'])
            self.slowDAQ_event = None

            for key, cb in self.slowDAQ_t0_checkbuttons.items():
                cb.config(state=DISABLED)
                self.slowDAQ_t0_vars[key].set(False)

            self.slowDAQ_ax.clear()
            self.slowDAQ_canvas.draw_idle()
            return

        path = os.path.join(self.raw_directory, self.run)

        try:
            self.slowDAQ_event = GetEvent(path, self.event, "run_control", "slow_daq")
        except FileNotFoundError:
            self.slowDAQ_error("No data")
            return
        except Exception as e:
            self.slowDAQ_error("GetEvent error", e)
            return

        try:
            data = self.slowDAQ_event.get('slow_daq', self.slowDAQ_event)

            sensor_keys = [
                k for k, v in data.items()
                if isinstance(v, np.ndarray) and v.ndim == 1
                and k not in ('time_ms', 'valves', 'loaded')
            ]
            sensor_keys.sort()

            previous = self.slowDAQ_combobox.get()
            self.slowDAQ_combobox['values'] = sensor_keys
            if sensor_keys:
                if previous in sensor_keys:
                    self.slowDAQ_combobox.set(previous)
                else:
                    self.slowDAQ_combobox.set(sensor_keys[0])
                self.slowDAQ_combobox.state(['!disabled', 'readonly'])
            else:
                self.slowDAQ_combobox.set('')
                self.slowDAQ_combobox.state(['disabled'])

            self.draw_slowDAQ()
        except Exception as e:
            self.slowDAQ_error("EventDisplay error", e)

        gc.collect()

    def process_slowDAQ_cuts(self, time_ms, y):
        # Ensure numpy arrays
        t = np.asarray(time_ms, dtype=float)
        x = np.asarray(y, dtype=float)

        if t.size == 0 or x.size == 0:
            return t, x

        # Time window
        tmin = self.slowDAQ_tmin if self.slowDAQ_tmin is not None else t[0]
        tmax = self.slowDAQ_tmax if self.slowDAQ_tmax is not None else t[-1]

        time_mask = (t >= tmin) & (t <= tmax)
        if not np.any(time_mask):
            # If the window excludes everything, just fall back to full data
            time_mask = np.ones_like(t, dtype=bool)

        t = t[time_mask]
        x = x[time_mask]

        # Amplitude window
        vmin = self.slowDAQ_ymin
        vmax = self.slowDAQ_ymax

        # If neither is set, we’re done
        if vmin is None and vmax is None:
            return t, x

        amp_mask = np.ones_like(x, dtype=bool)
        if vmin is not None:
            amp_mask &= (x >= vmin)
        if vmax is not None:
            amp_mask &= (x <= vmax)

        if not np.any(amp_mask):
            # If amp cuts remove everything, just use time window
            return t, x

        t = t[amp_mask]
        x = x[amp_mask]

        return t, x


    def apply_slowDAQ_cuts(self):
        def parse(var):
            s = var.get().strip()
            return float(s) if s else None

        self.slowDAQ_tmin = parse(self.slowDAQ_tmin_var)
        self.slowDAQ_tmax = parse(self.slowDAQ_tmax_var)
        self.slowDAQ_ymin = parse(self.slowDAQ_ymin_var)
        self.slowDAQ_ymax  = parse(self.slowDAQ_ymax_var)

        self.draw_slowDAQ()


    def draw_slowDAQ(self, event=None):
        # Refresh checkbutton enable/disable first: the guards below return early
        # for events with no slowDAQ trace, which would otherwise leave every box
        # showing the previous event's availability.
        self.update_slowDAQ_t0_widgets()

        if self.slowDAQ_event is None:
            return

        data = self.slowDAQ_event.get('slow_daq', self.slowDAQ_event)

        time_ms = data.get('time_ms', None)
        if time_ms is None:
            self.slowDAQ_error("'time_ms' not found in slowDAQ event")
            return

        sensor_name = self.slowDAQ_combobox.get()
        if not sensor_name:
            return

        y = data.get(sensor_name, None)
        if y is None:
            self.slowDAQ_sensor_error(sensor_name)
            return

        n = min(len(time_ms), len(y))
        if n == 0:
            print(f"{sensor_name} and time_ms length mismatch")
            self.slowDAQ_sensor_error(sensor_name)
            return

        # Choose the time frame ONCE per draw so the trace, t0 lines and pulses
        # cannot disagree, then shift before cutting so the t min/t max entries
        # mean the same thing as the axis the user is reading.
        latch = self.slowDAQ_frame_latch()
        t_abs = np.asarray(time_ms[:n], dtype=float)
        if latch is not None:
            t_abs = t_abs - latch

        t_cut, y_cut = self.process_slowDAQ_cuts(t_abs, y[:n])

        # Plot
        self.slowDAQ_ax.clear()
        self.slowDAQ_ax.plot(t_cut, y_cut, label=sensor_name)
        self.slowDAQ_ax.set_xlabel('ms since trigger latch' if latch is not None
                                   else 'Time [ms]')
        self.slowDAQ_ax.set_title(f"{sensor_name} {self.run}-{self.event}")
        self.slowDAQ_ax.grid(True)

        self.draw_scint_lollipop(latch)
        self.draw_slowDAQ_t0_lines(latch)
        self.apply_slowDAQ_default_xlim(latch)

        self.slowDAQ_fig.tight_layout()
        self.slowDAQ_canvas.draw_idle()

    def slowDAQ_frame_latch(self):
        # Latch [ms] to shift the x axis by, or None to stay on absolute time_ms.
        # Without scint_t0.sbc there is no latch, so the trigger-relative frame is
        # simply unavailable and the tab falls back to the original absolute axis.
        try:
            ev = int(self.event)
        except (TypeError, ValueError):
            return None
        latch = self.scint_latch.get(ev)
        if latch is None or not np.isfinite(latch):
            return None
        return float(latch)

    def get_acous_window_ms(self):
        # (pre, post) acoustic window half-widths [ms] from the run's rc.json,
        # or None. This is the span the acoustic DAQ actually recorded, so it is
        # the natural default zoom rather than an arbitrary crop.
        try:
            rc_path = os.path.join(self.raw_directory, self.run, 'rc.json')
            with open(rc_path, mode='r', encoding='utf-8') as f:
                acous = json.load(f)['acous']
            return float(acous['pre_trig_len']) * 1000.0, float(acous['post_trig_len']) * 1000.0
        except Exception:
            return None

    def apply_slowDAQ_default_xlim(self, latch):
        # Default the view to the acoustic window. Only a view default: explicit
        # cuts win, and no data is dropped, so zooming out still shows everything.
        if latch is None:
            return
        if self.slowDAQ_tmin is not None or self.slowDAQ_tmax is not None:
            return
        window = self.get_acous_window_ms()
        if window is None:
            return
        pre, post = window
        self.slowDAQ_ax.set_xlim(-pre, post)

    def update_slowDAQ_t0_widgets(self):
        try:
            ev = int(self.event)
        except (TypeError, ValueError):
            ev = None

        for key, cb in self.slowDAQ_t0_checkbuttons.items():
            val = self.get_slowDAQ_t0(key, ev) if ev is not None else None
            if val is not None and np.isfinite(val):
                cb.config(state=NORMAL)
            else:
                cb.config(state=DISABLED)
                self.slowDAQ_t0_vars[key].set(False)

        # Mirror the t0 checkbuttons: disable when the event can't supply pulses,
        # rather than letting the user tick a box that silently unticks itself.
        if self.scint_lollipop_unavailable(ev) is None:
            self.slowDAQ_scint_checkbutton.config(state=NORMAL)
        else:
            self.slowDAQ_scint_checkbutton.config(state=DISABLED)
            self.slowDAQ_scint_var.set(False)

    def get_slowDAQ_t0(self, key, ev):
        if key == 'pressure':
            return self.pressure_t0.get(ev)
        elif key == 'scint':
            return self.scint_t0.get(ev)
        elif key == 'bubble':
            return self.bubble_t0_ms(ev)

    def draw_slowDAQ_t0_lines(self, latch=None):
        # Draw a dashed vertical line per enabled t0 type at its event value+offset
        try:
            ev = int(self.event)
        except (TypeError, ValueError):
            return

        # The t0 values are already trigger-relative. When the axis has been
        # shifted to latch = 0 they need no offset at all; only the fallback
        # absolute axis needs t_compression to place them.
        offset_auto = 0.0 if latch is not None else self.t_compression.get(ev, 0.0)

        for key, color in SLOWDAQ_T0_COLORS.items():
            if not self.slowDAQ_t0_vars[key].get():
                continue
            val = self.get_slowDAQ_t0(key, ev)
            if val is None or not np.isfinite(val):
                continue
            self.slowDAQ_ax.axvline(
                x=val + offset_auto,
                linestyle='dashed', color=color, label=f'{key} t0',
            )

        # Merge the twin axis's handles in; a legend built from slowDAQ_ax alone
        # would silently omit the scint pulses drawn on the overlay.
        handles, labels = self.slowDAQ_ax.get_legend_handles_labels()
        if self.slowDAQ_scint_ax is not None and self.slowDAQ_scint_ax.get_visible():
            h2, l2 = self.slowDAQ_scint_ax.get_legend_handles_labels()
            handles += h2
            labels += l2

        if handles:
            self.slowDAQ_ax.legend(handles, labels)

    def load_scint_lollipop(self, ev):
        # Per-CAEN-trigger (time, summed pulse area) for one event, or None.
        # Time comes from the raw file's TriggerTimeTag and area from the recon
        # file; recon carries no timestamps and raw carries no pulse areas, so
        # both are needed. Recon row N is raw trigger N.
        if ev in self._scint_lolli_cache:
            return self._scint_lolli_cache[ev]

        result = None
        try:
            if self._scint_lolli_run != self.run or self._scint_lolli_all is None:
                self._scint_lolli_all = None
                self._scint_lolli_cache = {}  # event numbers repeat across runs
                self._scint_lolli_run = self.run
                path = self._find_recon('scintillation.sbc', self.run)
                if path is not None:
                    self._scint_lolli_all = Streamer(path).data

            recon = self._scint_lolli_all
            if recon is not None:
                if 'ev' in recon.dtype.names:
                    recon = recon[recon['ev'] == ev]
                if len(recon):
                    # NaN marks a channel with no hit, so nansum gives 0 for a
                    # trigger with nothing on any channel.
                    area = np.nansum(np.asarray(recon['hit_area']), axis=1)

                    # TriggerTimeTag is read eagerly even in lazy mode, so this
                    # never pulls waveforms off disk.
                    run_path = os.path.join(self.raw_directory, self.run)
                    event = GetEvent(run_path, ev, 'run_control', 'scintillation',
                                     strictMode=False, lazy_load_scintillation=True)
                    scint = event['scintillation']
                    if scint.get('loaded'):
                        tt = _unwrap_caen_timestamp(scint['TriggerTimeTag'])
                        t_ms = (tt - tt[0]) / SCINT_SAMPLE_RATE * 1000.0
                        n = min(len(t_ms), len(area))
                        if n:
                            result = (np.asarray(t_ms[:n], dtype=float),
                                      np.asarray(area[:n], dtype=float))
        except Exception as e:
            self.logger.error(
                'failed to load scint pulses for ev {}: {}'.format(ev, e))
            result = None

        self._scint_lolli_cache[ev] = result
        return result

    def scint_lollipop_unavailable(self, ev):
        # Reason the overlay can't be drawn, or None if it can. Cheap checks only
        # (dict lookups + a path probe) so this is safe to call on every redraw.
        if ev is None:
            return 'no event selected'
        try:
            # reco_directory is set from the config after this tab is constructed
            if self._find_recon('scintillation.sbc', self.run) is None:
                return 'no recon scintillation.sbc for this run'
        except AttributeError:
            return 'reco directory not configured yet'
        latch = self.scint_latch.get(ev)
        if latch is None or not np.isfinite(latch):
            return 'no scint t0 for this event (missing scint_t0.sbc, or Failed)'
        return None

    def draw_scint_lollipop(self, latch=None):
        # Overlay one stem per CAEN trigger on a twin y axis. Pulse area and the
        # slowDAQ sensor have unrelated units, so they must not share a scale.
        if self.slowDAQ_scint_ax is not None:
            self.slowDAQ_scint_ax.clear()
            self.slowDAQ_scint_ax.set_visible(False)

        if not self.slowDAQ_scint_var.get():
            return

        try:
            ev = int(self.event)
        except (TypeError, ValueError):
            return

        reason = self.scint_lollipop_unavailable(ev)
        if reason is None and self.load_scint_lollipop(ev) is None:
            reason = 'recon/raw scintillation could not be read'
        if reason is not None:
            # Say so on the plot; a silently empty overlay reads as a broken checkbox.
            self.slowDAQ_ax.set_title(
                '{}  |  scint pulses: {}'.format(
                    self.slowDAQ_ax.get_title(), reason), fontsize='small')
            self.logger.info('scint lollipop unavailable: {}'.format(reason))
            return

        t_caen, area = self.load_scint_lollipop(ev)
        if latch is None:
            latch = self.scint_latch[ev]

        # t_caen is measured from the first CAEN trigger, and slowDAQ time_ms
        # shares that origin (both start at event data-taking), so subtracting the
        # latch puts the pulses and the trace on one axis with the latch at 0.
        # No t_compression is involved.
        x = t_caen - latch

        # Honour the tab's time cuts. An event carries tens of thousands of
        # triggers over the full expansion, so without this the stems are a
        # solid block; t min/t max is how the user zooms to the window of
        # interest around t0.
        if self.slowDAQ_tmin is not None or self.slowDAQ_tmax is not None:
            lo = self.slowDAQ_tmin if self.slowDAQ_tmin is not None else -np.inf
            hi = self.slowDAQ_tmax if self.slowDAQ_tmax is not None else np.inf
            mask = (x >= lo) & (x <= hi)
            x, area = x[mask], area[mask]
            if not len(x):
                return

        if self.slowDAQ_scint_ax is None:
            self.slowDAQ_scint_ax = self.slowDAQ_ax.twinx()
        ax2 = self.slowDAQ_scint_ax
        ax2.set_visible(True)

        # Encode area as colour as well as height, so a big pulse is obvious
        # without reading the axis. Normalise on log10 to match the symlog y
        # scale; areas span ~4 decades and a linear norm would wash out all but
        # the largest few. clip keeps the exact zeros inside the colour range.
        with np.errstate(divide='ignore'):
            log_area = np.log10(np.clip(area, 1.0, None))
        vmax = log_area.max() if len(log_area) and log_area.max() > 0 else 1.0
        colors = matplotlib.colormaps['seismic'](log_area / vmax)

        lc = ax2.vlines(x, 0.0, area, linewidth=0.5, alpha=0.5)
        lc.set_color(colors)
        ax2.scatter(x, area, s=5, c=log_area, cmap='seismic', vmin=0.0, vmax=vmax)
        # scatter with an array c= yields no usable legend handle, so carry a
        # proxy: without it the merged legend silently omits the pulses.
        ax2.plot([], [], linestyle='none', marker='o', markersize=4,
                 color=matplotlib.colormaps['seismic'](0.85),
                 label='scint pulse area (colour = log area)')
        # symlog keeps the many near-zero triggers readable without dropping the
        # non-positive areas a plain log scale would discard.
        ax2.set_yscale('symlog')
        ax2.set_ylabel('Pulse area summed over SiPMs')


    def slowDAQ_error(self, label, error=None):
        if error is not None:
            print(f"{label}: {error}")

        self.slowDAQ_event = None
        self.slowDAQ_combobox['values'] = []
        self.slowDAQ_combobox.set('')
        self.slowDAQ_combobox.state(['disabled'])

        self.slowDAQ_ax.clear()
        self.slowDAQ_ax.text(
            0.5, 0.5, f"{label} for {self.run} - {self.event}",
            transform=self.slowDAQ_ax.transAxes,
            ha='center', va='center', fontsize=12, wrap=True
        )
        self.slowDAQ_ax.set_xlabel("Time [ms]")
        self.slowDAQ_ax.set_title(f"Slow DAQ - {self.run}-{self.event}")
        self.slowDAQ_canvas.draw_idle()

    def slowDAQ_sensor_error(self, sensor_name):
        # show error when sensor data is corrupted or malformed
        message = f"'{sensor_name}' not found for {self.run}-{self.event}."
        print(message)
        self.slowDAQ_ax.clear()
        self.slowDAQ_ax.text(
            0.5, 0.5, message,
            transform=self.slowDAQ_ax.transAxes,
            ha='center', va='center', fontsize=12, wrap=True
        )
        self.slowDAQ_ax.set_xlabel("Time [ms]")
        self.slowDAQ_ax.set_title(f"{sensor_name} {self.run}-{self.event}")
        self.slowDAQ_canvas.draw_idle()

