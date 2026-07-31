# Imports
import gc
import os
import matplotlib
import scipy.signal
import tkinter as tk
from tkinter import ttk, DISABLED, NORMAL
import numpy as np
import sys

matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from PIL import Image, ImageTk
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from GetEvent import GetEvent
from t0_common import T0_COLORS, SIPM_AREA_PER_PHD

# Default zoom pad either side of the (pressure t0, trigger) pair, in seconds.
PIEZO_T0_ZOOM_PAD_S = 0.100

# Fraction of one trace-to-trace step reserved below the bottom trace for the
# SiPM pulse stems.
PIEZO_PULSE_BAND_FRAC = 0.8


class Piezo(tk.Frame):
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)

        # For the fastDAQ tab
        self.piezo_cutoff_low = 2000
        self.piezo_cutoff_high = 10000
        self.piezo_beginning_time = -.1
        self.piezo_ending_time = 0.0
        self.piezo_max_points = 4000
        self.piezo_timerange_checkbutton_var = tk.BooleanVar(value=True)

        # selected_piezos persists channel-name selections across events/runs.
        # Tries to match selections from the previous event, falls back to the
        # first available channel to avoid an empty plot.
        self.selected_piezos = None
        self.piezo_channels = []
        self.piezo_checkbox_vars = {}
        self.piezo_checkbox_widgets = []

        # Trigger-relative t0 overlays, all from the recon pipeline.
        self.piezo_t0_vars = {
            key: tk.BooleanVar(value=False) for key in T0_COLORS
        }
        self.piezo_t0_checkbuttons = {}
        self.piezo_scint_var = tk.BooleanVar(value=False)

        # Initial Functions
        self.create_piezo_widgets()
        self.piezo_canvas_setup()

    def piezo_canvas_setup(self):
        # These used to be created every time a new event was loaded or something
        # changed in the log viewer. This lead to many figs, ax, and canvas that
        # continued to fill memory until python crashed.

        # Create Figures, Axes, and Canvases
        self.piezo_fig = Figure(figsize=(7, 5), dpi=100)
        self.piezo_ax = self.piezo_fig.add_subplot(111)
        self.piezo_canvas = FigureCanvasTkAgg(self.piezo_fig, self.piezo_tab_right)

        # Matplotlib's built-in pan and zoom toolbar
        self.piezo_toolbar = NavigationToolbar2Tk(
            self.piezo_canvas, self.piezo_tab_right, pack_toolbar=False)
        self.piezo_toolbar.update()
        self.piezo_toolbar.grid(row=1, column=1, sticky='w')

    def load_fastDAQ_piezo(self):
        if not self.load_fastDAQ_piezo_checkbutton_var.get():
            self.piezo_tab_right.grid_forget()
            self.disable_piezo_t0_widgets()
            return
        else:
            self.piezo_tab_right.grid(row=0, column=1, sticky='NW')

        if self.zip_flag:
            path = os.path.join(self.raw_directory, self.run, '.zip')

        path = os.path.join(self.raw_directory, self.run)

        try:
            self.fastDAQ_event = GetEvent(path, self.event, "run_control", "acoustics",
                                          physical_units=True)
        except FileNotFoundError:
            self.piezo_error("No data")
            return
        except Exception as e:
            self.piezo_error("GetEvent error", e)
            return

        try:
            acous_config = self.fastDAQ_event['run_control'].get('acous', {})
            num_channels = self.fastDAQ_event['acoustics']['Waveforms_V'].shape[1]
            channels = [
                acous_config.get(f'ch{i+1}', {}).get('name', f'Channel {i+1}')
                for i in range(num_channels)
            ]
            self.piezo_channels = channels

            self._rebuild_piezo_checkboxes(channels)
            self.draw_fastDAQ_piezo()
        except Exception as e:
            self.piezo_error("EventDisplay error", e)

        # Garbage Collecting
        gc.collect()

    def _rebuild_piezo_checkboxes(self, channels):
        for w in self.piezo_checkbox_widgets:
            w.destroy()
        self.piezo_checkbox_widgets = []
        self.piezo_checkbox_vars = {}

        # Use previous selected channels if possible
        remembered = set(self.selected_piezos or [])
        for name in channels:
            var = tk.BooleanVar(value=(name in remembered))
            cb = tk.Checkbutton(
                self.piezo_channel_frame,
                text=name,
                variable=var,
                anchor='w')
            cb.pack(anchor='w', fill='x')
            self.piezo_checkbox_vars[name] = var
            self.piezo_checkbox_widgets.append(cb)

        # If nothing carried over default to the first channel.
        if channels and not any(var.get() for var in self.piezo_checkbox_vars.values()):
            self.piezo_checkbox_vars[channels[0]].set(True)

    def update_piezo_t0_widgets(self):
        # Enable only the overlays this event can actually supply, and clear the
        # var when disabling so a box never stays ticked with nothing behind it.
        try:
            ev = int(self.event)
        except (TypeError, ValueError):
            ev = None

        for key, cb in self.piezo_t0_checkbuttons.items():
            val = self.get_t0_ms(key, ev) if ev is not None else None
            if val is not None and np.isfinite(val):
                cb.config(state=NORMAL)
            else:
                cb.config(state=DISABLED)
                self.piezo_t0_vars[key].set(False)

        if ev is not None and self.scint_pulses_unavailable(ev) is None:
            self.piezo_scint_checkbutton.config(state=NORMAL)
        else:
            self.piezo_scint_checkbutton.config(state=DISABLED)
            self.piezo_scint_var.set(False)

    def disable_piezo_t0_widgets(self):
        for key, cb in self.piezo_t0_checkbuttons.items():
            cb.config(state=DISABLED)
            self.piezo_t0_vars[key].set(False)
        self.piezo_scint_checkbutton.config(state=DISABLED)
        self.piezo_scint_var.set(False)

    def draw_fastDAQ_piezo(self):
        # Refresh availability before the early returns below, or an event with no
        # piezo data leaves the boxes advertising the previous event's data.
        self.update_piezo_t0_widgets()

        if not self.load_fastDAQ_piezo_checkbutton_var.get():
            self.piezo_tab_right.grid_forget()
            self.disable_piezo_t0_widgets()
            return
        else:
            self.piezo_tab_right.grid(row=0, column=1, sticky='NW')

        self.piezo_cutoff_low = int(self.piezo_cutoff_low_entry.get())
        if(self.piezo_cutoff_low < 1):
            self.piezo_cutoff_low = 1
            self.piezo_cutoff_low_entry.delete(0,tk.END)
            self.piezo_cutoff_low_entry.insert(0,self.piezo_cutoff_low)
        self.piezo_cutoff_high = int(self.piezo_cutoff_high_entry.get())
        self.piezo_beginning_time = float(self.piezo_beginning_time_entry.get())
        self.piezo_ending_time = float(self.piezo_ending_time_entry.get())
        try:
            self.piezo_max_points = max(1, int(self.piezo_max_points_entry.get()))
        except ValueError:
            self.piezo_max_points = 4000
        self.piezo_max_points_entry.delete(0, tk.END)
        self.piezo_max_points_entry.insert(0, self.piezo_max_points)
        # Snapshot the current selection by name. If the user unchecked the
        # last one, snap back to channels[0]
        selected = [name for name, var in self.piezo_checkbox_vars.items() if var.get()]
        if not selected and self.piezo_channels:
            first = self.piezo_channels[0]
            self.piezo_checkbox_vars[first].set(True)
            selected = [first]
        self.selected_piezos = selected

        self.draw_filtered_piezo_trace(selected)

    def draw_filtered_piezo_trace(self, selected_names):
        # A failed load leaves fastDAQ_event None with its message already on the
        # axes. Any redraw after that (a checkbutton, 'reload') would subscript None,
        # which the handler below does not catch, so stop here and keep the message.
        if self.fastDAQ_event is None:
            return

        try:
            acoustics = self.fastDAQ_event['acoustics']
            piezo_time = np.asarray(acoustics['time_s'])

            # Nyquist frequency
            # scipy.signal.butter requires 0 < Wn < 1 strictly.
            # no need for lower clamp because we force the widget to be at least 1 hz
            # low-pass at cutoff_high, high-pass at cutoff_low
            fn = acoustics['sample_rate'] / 2
            max_wn = 0.99999
            high_wn = self.piezo_cutoff_high / fn
            if high_wn >= 1:
                self.logger.error('Cutoff high >= Nyquist, clamping below Nyquist')
                high_wn = max_wn
                self.piezo_cutoff_high = int(fn * high_wn)
                self.piezo_cutoff_high_entry.delete(0, tk.END)
                self.piezo_cutoff_high_entry.insert(0, self.piezo_cutoff_high)

            low_wn = self.piezo_cutoff_low / fn
            if low_wn >= high_wn or low_wn >= 1:
                self.logger.error('Cutoff low >= cutoff high, clamping')
                low_wn = high_wn / 2
                self.piezo_cutoff_low = max(int(fn * low_wn), 1)
                self.piezo_cutoff_low_entry.delete(0, tk.END)
                self.piezo_cutoff_low_entry.insert(0, self.piezo_cutoff_low)

            # Set Plot Labels
            self.piezo_ax.clear()
            self.piezo_ax.set_title(str(self.run) + " " + str(self.event))
            self.piezo_ax.set_xlabel('[s]')
            self.piezo_ax.set_yticks([])

            if not self.piezo_timerange_checkbutton_var.get():
                self.piezo_ending_time_entry['state'] = tk.NORMAL
                self.piezo_beginning_time_entry['state'] = tk.NORMAL
                self.piezo_beginning_time_label['state'] = tk.NORMAL
                self.piezo_ending_time_label['state'] = tk.NORMAL
                window = (piezo_time > self.piezo_beginning_time) & (piezo_time < self.piezo_ending_time)
                plot_time = piezo_time[window]
                self.piezo_ax.set_xlim(self.piezo_beginning_time, self.piezo_ending_time)
            else:
                self.piezo_ending_time_entry['state'] = tk.DISABLED
                self.piezo_beginning_time_entry['state'] = tk.DISABLED
                self.piezo_beginning_time_label['state'] = tk.DISABLED
                self.piezo_ending_time_label['state'] = tk.DISABLED
                window = None
                plot_time = piezo_time
                self.piezo_ax.set_xlim(*self.piezo_default_xlim(piezo_time))

            plot_time_ds = self._decimate(plot_time,
                                          max_points=self.piezo_max_points)

            # Stack selected channels vertically in checkbox/channel order.
            cmap = matplotlib.colormaps['tab10']
            traces = []
            for idx, name in enumerate(self.piezo_channels):
                if name not in selected_names:
                    continue

                raw = np.asarray(acoustics['Waveforms_V'][0][idx])
                filtered = self._bandpass_piezo(raw, low_wn, high_wn)

                if window is not None:
                    filtered = filtered[window]

                v_ds = self._decimate(filtered, max_points=self.piezo_max_points)
                if v_ds.size:
                    traces.append((name, cmap(idx % cmap.N), v_ds - np.mean(v_ds)))

            step = 1.2 * max((np.ptp(v) for _, _, v in traces), default=0.0)
            
            for i, (name, color, v) in enumerate(reversed(traces)):
                offset = i * step
                self.piezo_ax.plot(plot_time_ds, v + offset, color=color)
                # Annotate axis instead of legend
                self.piezo_ax.annotate(
                    name, xy=(0, offset), xycoords=('axes fraction', 'data'),
                    xytext=(-4, 0), textcoords='offset points',
                    color=color, va='center', ha='right', fontsize=8,
                    annotation_clip=False)

            self.piezo_ax.relim()
            self.piezo_ax.autoscale_view(scalex=False, scaley=True)

            # Overlays go after autoscale: the pulse band is sized from the trace
            # stack and then extends the y limits itself.
            self.draw_piezo_t0_lines()
            self.draw_piezo_scint_pulses(traces, step)

            # One legend for every overlay. The traces themselves are labelled by
            # the axis annotations above, so they contribute no handles.
            handles, labels = self.piezo_ax.get_legend_handles_labels()
            if handles:
                self.piezo_ax.legend(handles, labels, fontsize='small', loc='upper right')

            # Update Canvas
            self.piezo_canvas.draw_idle()
            self.piezo_canvas.get_tk_widget().grid(row=0, column=1)

        ##added same handling for IndexError for when the piezo is not in the given multiboard. May lose specificity
        except (KeyError, IndexError, AttributeError):
            self.error += 'piezo data not found\n'
            self.destroy_children(self.piezo_tab_right)
            canvas = tk.Canvas(self.piezo_tab_right, width=self.init_image_width, height=self.init_image_height)
            self.reset_zoom(canvas)

            ### draw not found image
            image = Image.open('notfound.jpeg')
            self.native_image_width, self.native_image_height = image.size
            image = image.resize((int(canvas.image_width), int(canvas.image_height)),
                                 self.antialias_checkbutton_var.get())
            image = image.crop((canvas.crop_left, canvas.crop_bottom, canvas.crop_right, canvas.crop_top))

            canvas.image = canvas.create_image(0, 0, anchor=tk.NW, image=None)
            canvas.photo = ImageTk.PhotoImage(image)
            canvas.itemconfig(canvas.image, image=canvas.photo)
            canvas.grid(row=0, column=1, sticky='NW')

    def piezo_event_number(self):
        try:
            return int(self.event)
        except (TypeError, ValueError):
            return None

    def piezo_default_xlim(self, piezo_time):
        # Default view for the full-time-window mode: frame the interval between the
        # pressure t0 and the trigger, padded either side, since that is where the
        # acoustic signal of interest sits. Falls back to the whole recorded trace
        # (the acoustic window) when there is no usable pressure t0.
        #
        # View only -- nothing is masked in this mode, so zooming out still shows
        # the full trace.
        full = (float(piezo_time[0]), float(piezo_time[-1]))

        ev = self.piezo_event_number()
        if ev is None:
            return full
        pt0 = self.get_t0_ms('pressure', ev)
        if pt0 is None or not np.isfinite(pt0):
            return full

        pt0_s = pt0 / 1000.0
        # The trigger is t = 0 on this axis; only zoom when pt0 precedes it, else
        # the window would be inverted or degenerate.
        if pt0_s >= 0.0:
            return full

        lo = max(full[0], pt0_s - PIEZO_T0_ZOOM_PAD_S)
        hi = min(full[1], PIEZO_T0_ZOOM_PAD_S)
        if not lo < hi:
            return full
        return (lo, hi)

    def draw_piezo_t0_lines(self):
        # Dashed vertical line per enabled t0. All three are milliseconds relative
        # to the trigger, which is exactly this axis's origin (see t0_common), so a
        # unit conversion is the only transform needed.
        ev = self.piezo_event_number()
        if ev is None:
            return

        for key, color in T0_COLORS.items():
            if not self.piezo_t0_vars[key].get():
                continue
            val = self.get_t0_ms(key, ev)
            if val is None or not np.isfinite(val):
                continue
            self.piezo_ax.axvline(x=val / 1000.0, linestyle='dashed', color=color,
                                  label=f'{key} t0')

    def draw_piezo_scint_pulses(self, traces, step):
        # One stem per CAEN trigger in a band below the bottom trace. Drawn on
        # piezo_ax rather than a twin axis: this axis has no y ticks and stacks
        # traces at arbitrary offsets, so a second scaled axis would have nothing
        # to line up against. Magnitude is carried by stem height plus colour, with
        # the absolute scale stated in the legend.
        #
        # A single row of stems, not one per trace, keeps a busy plot readable.
        if not self.piezo_scint_var.get():
            return

        ev = self.piezo_event_number()
        if ev is None:
            return

        reason = self.scint_pulses_unavailable(ev)
        pulses = None
        if reason is None:
            pulses = self.load_scint_pulses(ev)
            if pulses is None:
                reason = 'recon/raw scintillation could not be read'
        if reason is None and (not traces or step <= 0):
            reason = 'no piezo trace to anchor to'
        if reason is not None:
            # Say so on the plot; a silently empty overlay reads as a broken checkbox.
            self.piezo_ax.set_title(
                '{}  |  SiPM pulses: {}'.format(self.piezo_ax.get_title(), reason),
                fontsize='small')
            self.logger.info('piezo SiPM pulses unavailable: {}'.format(reason))
            return

        t_caen, area = pulses
        latch = self.scint_latch[ev]
        # t_caen is measured from the first CAEN trigger and latch is the acoustic
        # trigger in that same frame, so subtracting it lands on this axis directly.
        x = (t_caen - latch) / 1000.0

        # An event carries tens of thousands of triggers across the whole
        # expansion; without clipping to the view the band is a solid block.
        lo, hi = self.piezo_ax.get_xlim()
        mask = (x >= lo) & (x <= hi)
        x, area = x[mask], area[mask]
        if not len(x):
            return

        # Rescale to an approximate photo-electron count so stem heights sit on a
        # meaningful order of magnitude. See SIPM_AREA_PER_PHD: this is a ballpark,
        # not a calibration, hence the '~' on the legend below.
        phd = area / SIPM_AREA_PER_PHD

        # Areas span several decades, so scale on log10; a linear scale would wash
        # out everything but the largest few pulses. clip keeps zeros in range.
        log_phd = np.log10(np.clip(phd, 1.0, None))
        vmax = log_phd.max() if log_phd.max() > 0 else 1.0

        # traces is drawn reversed with offset i*step, so traces[-1] is the bottom
        # of the stack. Which channel that is changes as the user toggles channels,
        # hence recomputing it here instead of pinning a channel name.
        y0 = float(np.min(traces[-1][2]))
        band = PIEZO_PULSE_BAND_FRAC * step
        base = y0 - band
        tops = base + band * (log_phd / vmax)

        cmap = matplotlib.colormaps['seismic']
        lc = self.piezo_ax.vlines(x, base, tops, linewidth=0.5, alpha=0.5)
        lc.set_color(cmap(log_phd / vmax))
        self.piezo_ax.scatter(x, tops, s=5, c=log_phd, cmap='seismic',
                              vmin=0.0, vmax=vmax)
        # scatter with an array c= yields no usable legend handle, so carry a proxy;
        # it also holds the absolute scale the unlabelled band cannot show.
        self.piezo_ax.plot(
            [], [], linestyle='none', marker='o', markersize=4, color=cmap(0.85),
            label='SiPM pulses (colour/height = log area, max ~{:.2e} phd)'.format(phd.max()))

        # Make room for the band; autoscale has already run by this point.
        ymin, ymax = self.piezo_ax.get_ylim()
        self.piezo_ax.set_ylim(min(ymin, base - 0.1 * band), ymax)

    def _decimate(self, *arrays, max_points):
        n = len(arrays[0])
        stride = max(1, n // max_points)
        if stride == 1:
            return arrays if len(arrays) > 1 else arrays[0]
        out = tuple(a[::stride] for a in arrays)
        return out if len(out) > 1 else out[0]

    def _bandpass_piezo(self, v, low_wn, high_wn):
        b, a = scipy.signal.butter(3, high_wn)
        v = scipy.signal.lfilter(b, a, v)
        b, a = scipy.signal.butter(3, low_wn, 'high')
        return scipy.signal.lfilter(b, a, v)

    def destroy_children(self, frame):
        try:
            for widget in frame.winfo_children():
                widget.destroy()
        except AttributeError:
            pass

    def create_piezo_widgets(self):
        self.piezo_tab = tk.Frame(self.notebook)
        self.notebook.add(self.piezo_tab, text='Piezo')

        # Piezos tab
        # First setup frames for piezos tab
        self.piezo_tab_left = tk.Frame(self.piezo_tab, bd=5, relief=tk.SUNKEN)
        self.piezo_tab_left.grid(row=0, column=0, sticky='NW')

        self.piezo_tab_right = tk.Frame(self.piezo_tab, bd=5, relief=tk.SUNKEN)
        self.piezo_tab_right.grid(row=0, column=1, sticky='NW')

        #         self.piezo_scrollbar = tk.Scrollbar(self.piezo_tab_right, orient = 'vertical')
        #         self.piezo_scrollbar.pack(side = 'left', fill = 'y')
        #         self.piezo_scrollbar.grid(row = 0, column = 0, sticky = tk.N + tk.S + tk.W + tk.E)

        # Now within the piezos frames setup stuff
        self.load_fastDAQ_piezo_checkbutton = tk.Checkbutton(
            self.piezo_tab_left,
            text='Load fastDAQ',
            variable=self.load_fastDAQ_piezo_checkbutton_var,
            command=self.load_fastDAQ_piezo)
        self.load_fastDAQ_piezo_checkbutton.grid(row=0, column=0, columnspan=2, sticky='WE')

        ttk.Separator(self.piezo_tab_left, orient='horizontal').grid(row=1, column=0, columnspan=2, sticky='WE', pady=4)

        self.piezo_label = tk.Label(self.piezo_tab_left, text='Piezos:')
        self.piezo_label.grid(row=2, column=0, columnspan=2, sticky='W')

        # Container for the per-channel checkbuttons. Populated by
        # _rebuild_piezo_checkboxes each time an event loads.
        self.piezo_channel_frame = tk.Frame(self.piezo_tab_left)
        self.piezo_channel_frame.grid(row=3, column=0, columnspan=2, sticky='WE')

        ttk.Separator(self.piezo_tab_left, orient='horizontal').grid(row=4, column=0, columnspan=2, sticky='WE', pady=4)

        self.piezo_cutoff_low_label = tk.Label(self.piezo_tab_left, text='Freq cutoff low:')
        self.piezo_cutoff_low_label.grid(row=5, column=0, sticky='WE')

        self.piezo_cutoff_low_entry = tk.Entry(self.piezo_tab_left, width=12)
        self.piezo_cutoff_low_entry.insert(0, self.piezo_cutoff_low)
        self.piezo_cutoff_low_entry.grid(row=5, column=1, sticky='WE')

        self.piezo_cutoff_high_label = tk.Label(self.piezo_tab_left, text='Freq cutoff high:')
        self.piezo_cutoff_high_label.grid(row=6, column=0, sticky='WE')

        self.piezo_cutoff_high_entry = tk.Entry(self.piezo_tab_left, width=12)
        self.piezo_cutoff_high_entry.insert(0, self.piezo_cutoff_high)
        self.piezo_cutoff_high_entry.grid(row=6, column=1, sticky='WE')

        ttk.Separator(self.piezo_tab_left, orient='horizontal').grid(row=7, column=0, columnspan=2, sticky='WE', pady=4)

        self.piezo_max_points_label = tk.Label(self.piezo_tab_left, text='Plot samples:')
        self.piezo_max_points_label.grid(row=8, column=0, sticky='WE')

        self.piezo_max_points_entry = tk.Entry(self.piezo_tab_left, width=12)
        self.piezo_max_points_entry.insert(0, self.piezo_max_points)
        self.piezo_max_points_entry.grid(row=8, column=1, sticky='WE')

        self.piezo_beginning_time_label = tk.Label(self.piezo_tab_left, text='Beginning Time:')
        self.piezo_beginning_time_label.grid(row=9, column=0, sticky='WE')

        self.piezo_beginning_time_entry = tk.Entry(self.piezo_tab_left, width=12)
        self.piezo_beginning_time_entry.insert(0, self.piezo_beginning_time)
        self.piezo_beginning_time_entry.grid(row=9, column=1, sticky='WE')

        self.piezo_ending_time_label = tk.Label(self.piezo_tab_left, text='Ending Time:')
        self.piezo_ending_time_label.grid(row=10, column=0, sticky='WE')

        self.piezo_ending_time_entry = tk.Entry(self.piezo_tab_left, width=12)
        self.piezo_ending_time_entry.insert(0, self.piezo_ending_time)
        self.piezo_ending_time_entry.grid(row=10, column=1, sticky='WE')

        self.piezo_timerange_checkbutton = tk.Checkbutton(
            self.piezo_tab_left, text='Full time window',
            variable=self.piezo_timerange_checkbutton_var,
            command=self.draw_fastDAQ_piezo)
        self.piezo_timerange_checkbutton.grid(row=11, column=0, columnspan=2, sticky='WE')

        row = 12
        for key in T0_COLORS:
            cb = tk.Checkbutton(
                self.piezo_tab_left,
                text=f'Show {key} t0',
                variable=self.piezo_t0_vars[key],
                command=self.draw_fastDAQ_piezo,
                state=DISABLED)
            cb.grid(row=row, column=0, columnspan=2, sticky='WE')
            self.piezo_t0_checkbuttons[key] = cb
            row += 1

        self.piezo_scint_checkbutton = tk.Checkbutton(
            self.piezo_tab_left,
            text='Show SiPM pulses',
            variable=self.piezo_scint_var,
            command=self.draw_fastDAQ_piezo,
            state=DISABLED)
        self.piezo_scint_checkbutton.grid(row=row, column=0, columnspan=2, sticky='WE')
        row += 1

        self.reload_fastDAQ_piezo_button = tk.Button(self.piezo_tab_left, text='reload',
                                                     command=self.draw_fastDAQ_piezo)
        self.reload_fastDAQ_piezo_button.grid(row=row, column=0, columnspan=2, sticky='WE')

    def piezo_error(self, label, error=None):
        if error is not None:
            print(f"{label}: {error}")

        self.fastDAQ_event = None
        for w in self.piezo_checkbox_widgets:
            w.destroy()
        self.piezo_checkbox_widgets = []
        self.piezo_checkbox_vars = {}
        self.piezo_channels = []
        self.disable_piezo_t0_widgets()
        self.piezo_ax.clear()
        self.piezo_ax.text(0.2, 0.5, f"{label} for {self.run} - {self.event}", transform=self.piezo_ax.transAxes, fontsize=15)

        self.piezo_ax.set_xlabel('[s]')
        self.piezo_ax.set_ylabel('[V]')

        # Make sure elements are on canvas before calling draw_idle
        self.piezo_tab_right.grid(row=0, column=1, sticky='NW')
        self.piezo_canvas.get_tk_widget().grid(row=0, column=1, sticky='NW')
        self.piezo_canvas.draw_idle()
