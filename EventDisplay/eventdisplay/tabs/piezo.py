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
# matplotlib.use('TkAgg')
matplotlib.use('Agg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from PIL import Image, ImageTk
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
from GetEvent import GetEvent


class Piezo(tk.Frame):
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)

        # For the fastDAQ tab
        self.piezo_cutoff_low = 2000
        self.piezo_cutoff_high = 10000
        self.piezo_beginning_time = -.1
        self.piezo_ending_time = 0.0
        self.incremented_piezo_event = False
        self.piezo_timerange_checkbutton_var = tk.BooleanVar(value=False)
        self.t0 = None

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

    def load_fastDAQ_piezo(self):
        if not self.load_fastDAQ_piezo_checkbutton_var.get():
            self.piezo_tab_right.grid_forget()
            return
        else:
            self.piezo_tab_right.grid(row=0, column=1, sticky='NW')

        if self.zip_flag:
            path = os.path.join(self.raw_directory, self.run, '.zip')

        path = os.path.join(self.raw_directory, self.run)

        try:
            selected = ["run_control", "acoustics"]
            self.fastDAQ_event = GetEvent(path, self.event, *selected)
            
            wf_key = "Waveforms" if "Waveforms" in self.fastDAQ_event["acoustics"] else "Waveform"
            channels = [f"Channel {i+1}" for i in range(self.fastDAQ_event['acoustics'][wf_key].shape[1])]
            self.piezo_combobox['values'] = channels

            if channels:
                self.piezo_combobox.set(channels[0])
                self.piezo_combobox.state(['!disabled', 'readonly'])
            else:
                self.piezo_combobox.set('')
                self.piezo_combobox.state(['disabled'])

            self.draw_fastDAQ_piezo()
        except Exception as e:
            print(e)
            self.piezo_error()

        # Garbage Collecting
        gc.collect()

    def check_t0_exists(self):
        try:
            self.t0 = self.reco_row['fastDAQ_t0']
            self.piezo_plot_t0_checkbutton.config(state=NORMAL)
        except:
            self.piezo_plot_t0_checkbutton.config(state=DISABLED)
            self.t0 = None

    def draw_fastDAQ_piezo(self):
        if not self.load_fastDAQ_piezo_checkbutton_var.get():
            self.piezo_tab_right.grid_forget()
            return
        else:
            self.piezo_tab_right.grid(row=0, column=1, sticky='NW')

        # if int(self.run_type) == 10:
        #     # logger.error('not allowed for run_type=10')
        #     self.error += 'not allowed to view piezo data for run_type=10\n'
        #     self.destroy_children(self.piezo_tab_right)
        #     return
        self.check_t0_exists()
        self.piezo = self.piezo_combobox.get()
        self.piezo_cutoff_low = int(self.piezo_cutoff_low_entry.get())
        if(self.piezo_cutoff_low < 1):
            self.piezo_cutoff_low = 1
            self.piezo_cutoff_low_entry.delete(0,tk.END)
            self.piezo_cutoff_low_entry.insert(0,self.piezo_cutoff_low)
        self.piezo_cutoff_high = int(self.piezo_cutoff_high_entry.get())
        self.piezo_beginning_time = float(self.piezo_beginning_time_entry.get())
        self.piezo_ending_time = float(self.piezo_ending_time_entry.get())
        self.draw_filtered_piezo_trace(self.piezo)

    def draw_filtered_piezo_trace(self, piezo):
        try:
            wf_key = "Waveforms" if "Waveforms" in self.fastDAQ_event["acoustics"] else "Waveform"
            piezo_v = self.fastDAQ_event['acoustics'][wf_key][0][self.piezo_combobox.current()]
            piezo_time = np.arange(len(piezo_v)) * (1 / self.fastDAQ_event['acoustics']['sample_rate'])
            fn = len(piezo_v)/(piezo_time[-1]-piezo_time[0])/2
            # if (self.piezo_cutoff_high / fn) > 1:
            #     self.logger.error('Cutoff freq > Nyquist, setting = Nyquist')
            #     self.piezo_cutoff_high = int(fn)
            #     self.piezo_cutoff_high_entry.delete(0,tk.END)
            #     self.piezo_cutoff_high_entry.insert(0,self.piezo_cutoff_high)
                
            # b, a = scipy.signal.butter(3, self.piezo_cutoff_high / fn)
            # filtered_piezo_v = scipy.signal.lfilter(b, a, piezo_v)
            # b, a = scipy.signal.butter(3, self.piezo_cutoff_low / fn, 'high')
            # filtered_piezo_v = scipy.signal.lfilter(b, a, filtered_piezo_v)

            filtered_piezo_v = piezo_v

            # Set Plot Labels
            self.piezo_ax.clear()
            self.piezo_ax.set_title(piezo + " " + str(self.run) + " " + str(self.event))
            self.piezo_ax.set_xlabel('[s]')
            self.piezo_ax.set_ylabel('[V]')

            # if not self.piezo_timerange_checkbutton_var.get():
            #     self.piezo_ending_time_entry['state'] = tk.NORMAL
            #     self.piezo_beginning_time_entry['state'] = tk.NORMAL
            #     self.piezo_beginning_time_label['state'] = tk.NORMAL
            #     self.piezo_ending_time_label['state'] = tk.NORMAL
            #     window = (piezo_time > self.piezo_beginning_time) & (piezo_time < self.piezo_ending_time)
            #     piezo_time = piezo_time[window]
            #     filtered_piezo_v = filtered_piezo_v[window]
            #     self.piezo_ax.set_xlim(self.piezo_beginning_time, self.piezo_ending_time)
            # else:
            #     self.piezo_ending_time_entry['state'] = tk.DISABLED
            #     self.piezo_beginning_time_entry['state'] = tk.DISABLED
            #     self.piezo_beginning_time_label['state'] = tk.DISABLED
            #     self.piezo_ending_time_label['state'] = tk.DISABLED
            #     self.piezo_ax.set_xlim(piezo_time[0], piezo_time[-1])

            # Plot - Create Lines First Time, Update Data For Line Afterwards
            try:
                self.piezo_ax.lines[0].set_xdata(piezo_time)
                self.piezo_ax.lines[0].set_ydata(filtered_piezo_v)
            except:
                # self.piezo_ax.plot(piezo_time, filtered_piezo_v)
                # We don't need to draw every data point unless with have a million-pixel display. shoot for 4000 data points max
                self.piezo_ax.plot(piezo_time[::int(len(piezo_time)/40)], filtered_piezo_v[::int(len(filtered_piezo_v)/40)])
            # Rescale Axis
            self.piezo_ax.relim()
            self.piezo_ax.autoscale_view()

            # Add line at t0
            #self.check_t0_exists()
            if self.piezo_plot_t0_checkbutton_var.get():
                if self.reco_row and self.t0:
                    self.piezo_ax.axvline(x=self.t0, linestyle='dashed', color='r', label='t0')
                    self.incremented_piezo_event = True
                else:
                    if self.incremented_piezo_event:
                        self.error += 't0 unavailable: no reco data found for current event.'
                    else:
                        self.logger.error('t0 unavailable: no reco data found for current event.')
                    self.piezo_plot_t0_checkbutton_var.set(False)
                    self.incremented_piezo_event = False

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

        self.piezo_label = tk.Label(self.piezo_tab_left, text='Piezo:')
        self.piezo_label.grid(row=1, column=0, sticky='WE')

        self.piezo_combobox = ttk.Combobox(self.piezo_tab_left, width=12)
        self.piezo_combobox.grid(row=1, column=1, sticky='WE')

        self.piezo_cutoff_low_label = tk.Label(self.piezo_tab_left, text='Freq cutoff low:')
        self.piezo_cutoff_low_label.grid(row=2, column=0, sticky='WE')

        self.piezo_cutoff_low_entry = tk.Entry(self.piezo_tab_left, width=12)
        self.piezo_cutoff_low_entry.insert(0, self.piezo_cutoff_low)
        self.piezo_cutoff_low_entry.grid(row=2, column=1, sticky='WE')

        self.piezo_cutoff_high_label = tk.Label(self.piezo_tab_left, text='Freq cutoff high:')
        self.piezo_cutoff_high_label.grid(row=3, column=0, sticky='WE')

        self.piezo_cutoff_high_entry = tk.Entry(self.piezo_tab_left, width=12)
        self.piezo_cutoff_high_entry.insert(0, self.piezo_cutoff_high)
        self.piezo_cutoff_high_entry.grid(row=3, column=1, sticky='WE')

        self.piezo_timerange_checkbutton = tk.Checkbutton(
            self.piezo_tab_left, text='Full time window',
            variable=self.piezo_timerange_checkbutton_var,
            command=self.draw_fastDAQ_piezo)
        self.piezo_timerange_checkbutton.grid(row=6, column=0, columnspan=2, sticky='WE')

        self.piezo_beginning_time_label = tk.Label(self.piezo_tab_left, text='Beginning Time:')
        self.piezo_beginning_time_label.grid(row=4, column=0, sticky='WE')

        self.piezo_beginning_time_entry = tk.Entry(self.piezo_tab_left, width=12)
        self.piezo_beginning_time_entry.insert(0, self.piezo_beginning_time)
        self.piezo_beginning_time_entry.grid(row=4, column=1, sticky='WE')

        self.piezo_ending_time_label = tk.Label(self.piezo_tab_left, text='Ending Time:')
        self.piezo_ending_time_label.grid(row=5, column=0, sticky='WE')

        self.piezo_ending_time_entry = tk.Entry(self.piezo_tab_left, width=12)
        self.piezo_ending_time_entry.insert(0, self.piezo_ending_time)
        self.piezo_ending_time_entry.grid(row=5, column=1, sticky='WE')

        self.piezo_plot_t0_checkbutton = tk.Checkbutton(
            self.piezo_tab_left,
            text='Show t0',
            variable=self.piezo_plot_t0_checkbutton_var,
            command=self.draw_fastDAQ_piezo)
        self.piezo_plot_t0_checkbutton.grid(row=7, column=0, columnspan=2, sticky='WE')

        self.reload_fastDAQ_piezo_button = tk.Button(self.piezo_tab_left, text='reload',
                                                     command=self.draw_fastDAQ_piezo)
        self.reload_fastDAQ_piezo_button.grid(row=8, column=0, sticky='WE')

    def piezo_error(self):
        print(f"No acoustics.sbc for {self.run} - {self.event}")
        self.fastDAQ_event = None
        self.piezo_combobox['values'] = []
        self.piezo_combobox.set('')
        self.piezo_ax.clear()
        self.piezo_ax.text(0.2, 0.5, f"No data for {self.run} - {self.event}", transform=self.piezo_ax.transAxes, fontsize=15)

        self.piezo_ax.set_xlabel('[s]')
        self.piezo_ax.set_ylabel('[V]')

        # Make sure elements are on canvas before calling draw_idle
        self.piezo_tab_right.grid(row=0, column=1, sticky='NW')
        self.piezo_canvas.get_tk_widget().grid(row=0, column=1, sticky='NW')
        self.piezo_canvas.draw_idle()
