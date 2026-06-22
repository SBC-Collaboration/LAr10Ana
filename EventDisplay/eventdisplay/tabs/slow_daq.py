# Imports
import gc
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


SLOWDAQ_T0_COLORS = {'pressure': 'r', 'bubble': 'g', 'scint': 'b'}


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

        # "Show t0" toggles. Disabled until a finite t0 value exists for the current event
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

        t_cut, y_cut = self.process_slowDAQ_cuts(time_ms[:n], y[:n])

        # Plot
        self.slowDAQ_ax.clear()
        self.slowDAQ_ax.plot(t_cut, y_cut)
        self.slowDAQ_ax.set_xlabel("Time [ms]")
        self.slowDAQ_ax.set_title(f"{sensor_name} {self.run}-{self.event}")
        self.slowDAQ_ax.grid(True)

        self.update_slowDAQ_t0_widgets()
        self.draw_slowDAQ_t0_lines()

        self.slowDAQ_fig.tight_layout()
        self.slowDAQ_canvas.draw_idle()

    def update_slowDAQ_t0_widgets(self):
        # Only enable each "show t0" checkbox only when a t0 exists for the current event
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

    def get_slowDAQ_t0(self, key, ev):
        if key == 'pressure':
            return self.pressure_t0.get(ev)
        elif key == 'scint':
            return self.scint_t0.get(ev)
        elif key == 'bubble':
            return self.bubble_t0_ms(ev)

    def draw_slowDAQ_t0_lines(self):
        # Draw a dashed vertical line per enabled t0 type at its event value+offset
        try:
            ev = int(self.event)
        except (TypeError, ValueError):
            return

        # t_compression maps the trigger t0 onto the absolute slowDAQ time_ms axis
        offset_auto = self.t_compression.get(ev, 0.0)

        drawn = False
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
            drawn = True

        if drawn:
            self.slowDAQ_ax.legend()
    
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

