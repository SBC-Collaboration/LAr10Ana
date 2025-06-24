# Imports
import gc
import os
import matplotlib
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
#
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from PICOcode.DataHandling.GetEvent import Event
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from PIL import Image, ImageTk
from scipy.optimize import curve_fit
import numpy as np
from scipy.signal import butter, sosfilt
from scipy import signal
import pandas as pd

class Dytran(tk.Frame):
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)

        # For the Dytran
        self.dytran_start_time = -0.10
        self.dytran_ending_time = 0.05
        self.dytran_curvefit_start_time = -0.05
        self.dytran_curvefit_window_time = 0.06
        self.dytran_default_start_time = self.dytran_start_time
        self.dytran_default_ending_time = self.dytran_ending_time
        self.dytran_curvefit_default_start_time = self.dytran_curvefit_start_time
        self.dytran_curvefit_default_window_time = self.dytran_curvefit_window_time
        self.incremented_dytran_event = False
        self.dytran_plot_t0_checkbutton_var = tk.BooleanVar(value=False)
        self.dytran_curvefit_checkbutton_var = tk.BooleanVar(value=False)
        self.dytran_timerange_checkbutton_var = tk.BooleanVar(value=False)
        self.dytran_lowpass_checkbutton_var = tk.BooleanVar(value=False)
        self.dytran_use_t0_curvefit_checkbutton_var = tk.BooleanVar(value=False)
        self.dytran_calc_t_end_curvefit_checkbutton_var = tk.BooleanVar(value=False)

        # Initial Functions
        self.create_dytran_widgets()
        self.dytran_canvas_setup()

    def dytran_canvas_setup(self):
        # These used to be created every time a new event was loaded or something
        # changed in the log viewer. This lead to many figs, ax, and canvas that
        # continued to fill memory until python crashed.

        # Create Figures, Axes, and Canvases
        self.dytran_fig = Figure(figsize=(7, 5), dpi=100)
        self.dytran_ax = self.dytran_fig.add_subplot(111)
        self.dytran_canvas = FigureCanvasTkAgg(self.dytran_fig, self.dytran_tab_right)

    def destroy_children(self, frame):
        try:
            for widget in frame.winfo_children():
                widget.destroy()
        except AttributeError:
            pass

    def load_fastDAQ_dytran(self):
        if not self.load_dytran_checkbutton_var.get():
            self.dytran_tab_right.grid_forget()
            return
        else:
            self.dytran_tab_right.grid(row=0, column=1, sticky='NW')

        if self.zip_flag:
            path = os.path.join(self.raw_directory, self.run, '.zip')
            
        path = os.path.join(self.raw_directory, self.run)
        self.fastDAQ_event = Event(path, self.event, ['fastDAQ'])
        
        self.dytran_combobox['values'] = [x for x in dir(self.fastDAQ_event.fastDAQ) if ('Dytran' in x) and ('_' not in x)]

        self.draw_fastDAQ_dytran()
        
        # Garbage Collecting
        gc.collect()

    def draw_fastDAQ_dytran(self):
        if not self.load_dytran_checkbutton_var.get():
            self.dytran_tab_right.grid_forget()
            return
        else:
            self.dytran_tab_right.grid(row=0, column=1, sticky='NW')

        if self.dytran_start_time_entry.get():
            self.dytran_start_time = float(self.dytran_start_time_entry.get())
        else:
            self.dytran_start_time = self.dytran_default_start_time
            messagebox.showwarning('Alert','Using default start time.' )
            self.dytran_start_time_entry.delete(0, tk.END)     
            self.dytran_start_time_entry.insert(0, self.dytran_start_time)

        if self.dytran_ending_time_entry.get():
            self.dytran_ending_time = float(self.dytran_ending_time_entry.get())
        else:
            self.dytran_ending_time = self.dytran_default_ending_time
            messagebox.showwarning('Alert','Using default end time.' )
            self.dytran_ending_time_entry.delete(0, tk.END)     
            self.dytran_ending_time_entry.insert(0, self.dytran_ending_time)

        if self.dytran_curvefit_start_time_entry.get():
            self.dytran_curvefit_start_time = float(self.dytran_curvefit_start_time_entry.get())
        else:
            self.dytran_curvefit_start_time = self.dytran_curvefit_default_start_time
            messagebox.showwarning('Alert','Using default fit start time.' )
            self.dytran_curvefit_start_time_entry.delete(0, tk.END)     
            self.dytran_curvefit_start_time_entry.insert(0, self.dytran_curvefit_start_time)

        if self.dytran_curvefit_window_time_entry.get():
            self.dytran_curvefit_window_time = float(self.dytran_curvefit_window_time_entry.get())
        else:
            self.dytran_curvefit_window_time = self.dytran_curvefit_default_window_time
            messagebox.showwarning('Alert','Using default fit window time.' )
            self.dytran_curvefit_window_time_entry.delete(0, tk.END)     
            self.dytran_curvefit_window_time_entry.insert(0, self.dytran_curvefit_window_time)

        self.dytran = self.dytran_combobox.get()
        self.draw_dytran_trace(self.dytran)
        
    def func(self, x, a, t0, c):
        return a*(x>=t0)*(x-t0)**4 + c

    def butter_lowpass_filter(self, data, lowcut, fs, order=6):
        nyq = 0.5 * fs
        low = lowcut / nyq
        if low <= 0:
            low = 0.0001
        sos = butter(order, low, analog=False, btype='lowpass', output='sos')
        y = sosfilt(sos, data)
        return y
    
    def draw_dytran_trace(self, dytran):
        board = 0
        if True:# try:
            #             if dytran not in self.fastDAQ_event['fastDAQ']['multiboards'][0]:
            #                 if dytran in self.fastDAQ_event['fastDAQ']['multiboards'][1]:
            #                     board = 1
            dytran_v = getattr(getattr(self.fastDAQ_event.fastDAQ, dytran), 'V')
            dytran_time = self.fastDAQ_event.fastDAQ.time[0]
            
            compression_std_def_hard = 10.0
            if self.dytran_lowpass_checkbutton_var.get():
                fs = int( 1.0 / (dytran_time[1] - dytran_time[0]) )
                dytran_v = self.butter_lowpass_filter(dytran_v, 1000, fs, 6)
                
            self.dytran_ax.clear()
            self.dytran_ax.set_title(dytran)
            self.dytran_ax.set_xlabel('[s]')
            self.dytran_ax.set_ylabel('[V]')
            if not self.dytran_timerange_checkbutton_var.get():
                self.dytran_ending_time_entry['state'] = tk.NORMAL
                self.dytran_start_time_entry['state'] = tk.NORMAL
                self.dytran_start_time_label['state'] = tk.NORMAL
                self.dytran_ending_time_label['state'] = tk.NORMAL
                # self.dytran_ax.set_xlim(self.dytran_start_time, self.dytran_ending_time)
                window = (dytran_time >= self.dytran_start_time) & (dytran_time <= self.dytran_ending_time)
            else:
                self.dytran_ending_time_entry['state'] = tk.DISABLED
                self.dytran_start_time_entry['state'] = tk.DISABLED
                self.dytran_start_time_label['state'] = tk.DISABLED
                self.dytran_ending_time_label['state'] = tk.DISABLED
                # self.dytran_ax.set_xlim(dytran_time[0], dytran_time[-1])
                window = (dytran_time > -999) & True

            # Plot - Create Lines First Time, Update Data For Line Afterwards
            try:
                self.dytran_ax.lines[0].set_xdata(dytran_time[window])
                self.dytran_ax.lines[0].set_ydata(dytran_v[window])
            except:
                # self.dytran_ax.plot(dytran_time[window], dytran_v[window])
                # We don't need to draw every data point unless with have a million-pixel display. shoot for 40000 data points max
                self.dytran_ax.plot(dytran_time[window][::int(len(dytran_time)/40000)], dytran_v[window][::int(len(dytran_v)/40000)])


            # Rescale Axis
            self.dytran_ax.relim()
            self.dytran_ax.autoscale_view()

            # Add line at t0
            if self.dytran_plot_t0_checkbutton_var.get():
                if self.reco_row:
                    self.dytran_ax.axvline(x=self.reco_row['fastDAQ_t0'], linestyle='dashed', color='r', label='t0')
                    self.incremented_dytran_event = True
                else:
                    if self.incremented_dytran_event:
                        self.error += 't0 unavailable: no reco data found for current event. Dytran trace not drawn'
                    else:
                        self.logger.error('t0 unavailable: no reco data found for current event. Dytran trace not drawn')
                    self.dytran_plot_t0_checkbutton_var.set(False)
                    self.incremented_dytran_event = False
                    
            if self.dytran_curvefit_checkbutton_var.get():
                self.dytran_curvefit_start_time_entry.config(state = tk.NORMAL)
                self.dytran_curvefit_window_time_entry.config(state = tk.NORMAL)
                self.dytran_calc_t_end_curvefit_checkbutton.config(state = tk.NORMAL)
                if self.reco_row:
                   self.dytran_use_t0_curvefit_checkbutton.config(state = tk.NORMAL)
                else:
                    self.dytran_curvefit_start_time = self.dytran_curvefit_default_start_time
                    self.dytran_curvefit_start_time_entry.delete(0, tk.END)     
                    self.dytran_curvefit_start_time_entry.insert(0, self.dytran_curvefit_start_time)
                    self.dytran_use_t0_curvefit_checkbutton.config(state = tk.DISABLED)
                
                if self.reco_row and self.dytran_use_t0_curvefit_checkbutton_var.get():
                    start_time = self.reco_row['fastDAQ_t0']
                    self.dytran_curvefit_start_time_entry.delete(0, tk.END)     
                    self.dytran_curvefit_start_time_entry.insert(0, start_time)
                    self.dytran_curvefit_start_time_entry.config(state = tk.DISABLED)
                else:
                   start_time = float(self.dytran_curvefit_start_time_entry.get())

                start_time_indx = np.where(dytran_time >= start_time)
                # print('start_time:', start_time)
                # print('start_time_indx:', start_time_indx)
                if len(start_time_indx[0]) == 0:         
                   lower_bound = 0
                else:
                   lower_bound = start_time_indx[0][0]                
                start_time = dytran_time[lower_bound]

                if self.dytran_calc_t_end_curvefit_checkbutton_var.get():
                    # use filtered dytran to do this
                    if not self.dytran_lowpass_checkbutton_var.get():
                        fs = int( 1.0 / (dytran_time[1] - dytran_time[0]) )
                        df = pd.Series( self.butter_lowpass_filter(dytran_v, 1000, fs, 6) )
                    else:
                        df = pd.Series(dytran_v)
                    dytran_std = df.rolling(window=1000).std()
                    if lower_bound == 0:
                        uncompressed_dytran_std = dytran_std[1000]
                    else:
                        uncompressed_dytran_std = np.mean(dytran_std[1000:lower_bound])
                    compression_indx = np.where(dytran_std > compression_std_def_hard * uncompressed_dytran_std)
                    # print('max, min, pre std: ', np.max(dytran_std), np.min(dytran_std), uncompressed_dytran_std)
                    if (len(compression_indx[0]) == 0):
                        messagebox.showerror('Error','Compression time could not be calculated.\nUsing previous fit window.')
                        end_time = start_time + float(self.dytran_curvefit_window_time_entry.get())
                    else:
                        end_time = dytran_time[compression_indx[0][0]] - 0.008
                        # print(' ', compression_indx[0][0], dytran_time[compression_indx[0][0]], dytran_std[compression_indx[0][0]])
                        self.dytran_curvefit_window_time_entry.delete(0, tk.END)     
                        self.dytran_curvefit_window_time_entry.insert(0, end_time - start_time)
                        self.dytran_curvefit_window_time_entry.config(state = tk.DISABLED)
                    # self.dytran_ax.plot(dytran_time[:], dytran_std, 'r')
                else:
                    end_time = start_time + float(self.dytran_curvefit_window_time_entry.get())

                end_time_indx = np.where(dytran_time >= end_time)                
                if len(end_time_indx[0]) == 0:     
                   upper_bound = -1 
                else:
                   upper_bound = end_time_indx[0][0]                 
                end_time = dytran_time[upper_bound]
                  
                if end_time < start_time:
                    messagebox.showerror('Error','The end time is lower than the start time. Setting default window.')
                    self.dytran_curvefit_window_time_entry.config(state = tk.NORMAL)
                    self.dytran_curvefit_window_time = self.dytran_curvefit_default_window_time
                    self.dytran_curvefit_window_time_entry.delete(0, tk.END)
                    self.dytran_curvefit_window_time_entry.insert(0, self.dytran_curvefit_window_time)
                    self.dytran_calc_t_end_curvefit_checkbutton_var.set(False)
                    
                try:
                    step = dytran_time[1] - dytran_time[0]
                    # if lower_bound < 100:
                    #     pre_t0_stddev = np.std(dytran_v[0:100])
                    # else:
                    #     pre_t0_stddev = np.std(dytran_v[0:lower_bound])
                    popt, pcov = curve_fit(self.func, dytran_time[lower_bound:upper_bound],  dytran_v[lower_bound:upper_bound], p0 = [1000, start_time, 0], bounds = ([0 , -0.20 , -10],[100000000, 0.20 , 10]))
                    # popt, pcov = curve_fit(self.func, dytran_time[lower_bound:upper_bound],  dytran_v[lower_bound:upper_bound], sigma=pre_t0_stddev*np.ones(len(dytran_v[lower_bound:upper_bound])), absolute_sigma=True, p0 = [1000, start_time, 0], bounds = ([0 , -0.20 , -10],[100000000, 0.20 , 10]))
                    a, t0 , c = popt
                    a_err, t0_err, c_err = np.sqrt(np.diag(pcov))
                    # determine goodness-of-fit R2
                    residuals = dytran_v[lower_bound:upper_bound] - self.func(dytran_time[lower_bound:upper_bound], *popt)
                    ss_res = np.sum(residuals**2)
                    ss_tot = np.sum((dytran_v[lower_bound:upper_bound]-np.mean(dytran_v[lower_bound:upper_bound]))**2)
                    r_squared = 1 - (ss_res / ss_tot)
                    text_comment = ''
                    if r_squared < 0.4:
                        text_comment = "###Poor fit###"
                    elif np.abs(t0) > 0.1:
                        text_comment = "###Weird t0###"
                    # print('info for run/ev:', self.run, self.event)
                    # print('r sq:', r_squared)
                    # print('popt here', popt)
                    # print('perr here', a_err, t0_err, c_err)
                    # print('pre-t0 stddev', pre_t0_stddev)
                    self.dytran_ax.text(-0.1 , 1.1 , "Curve fit to a(x - t0)^4 + c:" + f"\na = {a:.2E}\nt0 = {t0:.3f}\nc = {c:.2f}\n{text_comment}", horizontalalignment='left', verticalalignment='top' , color = '#b300b3', backgroundcolor = 'w', size = 15, transform = self.dytran_ax.transAxes)
                    # self.dytran_ax.text(0.04 , .85 , "Curve fit to a(x - t0)^4 + c :" + "\na = " +  f"{float(a):.2E} +- {float(a_err):.2E}" + '\nt0 = ' +  f"{float(t0):.2E} +- {float(t0_err):.2E}", horizontalalignment='left', verticalalignment='bottom' , color = '#b300b3', backgroundcolor = 'w', size = 15, transform = self.dytran_ax.transAxes)
                    # self.dytran_ax.axvline(x = start_time, linestyle='dotted', color='#b300b3', label='t_start')
                    # self.dytran_ax.axvline(x = end_time, linestyle='dotted', color='#b300b3', label='t_end')
                    time_cut = (dytran_time >= start_time) & (dytran_time <= end_time)
                    self.dytran_ax.plot(dytran_time[time_cut][::int(len(dytran_time)/40000)], dytran_v[time_cut][::int(len(dytran_v)/40000)], color='#b300b3')
                    x_line = np.arange(start_time , end_time , step)
                    y_line = self.func(x_line, a, t0, c)
                    # print('y = %.5f * (x-%3.3f)^4 + %3.3f' % (a, t0, c))
                    self.dytran_ax.plot(x_line, y_line, '--k')
                    
                except Exception as e:
                    messagebox.showerror('Error', 'Curvefit failed: ' + str(e) + '\nPerhaps change the fit window')

            else:
                self.dytran_curvefit_start_time_entry.config(state = tk.DISABLED)
                self.dytran_curvefit_window_time_entry.config(state = tk.DISABLED)                
                self.dytran_use_t0_curvefit_checkbutton.config(state = tk.DISABLED)
                self.dytran_calc_t_end_curvefit_checkbutton.config(state = tk.DISABLED)
                
            # Update Canvas
            self.dytran_canvas.draw()
            self.dytran_canvas.get_tk_widget().grid(row=0, column=1)
        else:#except (KeyError, IndexError):
            self.destroy_children(self.dytran_tab_right)
            self.error += 'dytran data not found\n'
            self.destroy_children(self.dytran_tab_right)
            canvas = tk.Canvas(self.dytran_tab_right, width=self.init_image_width, height=self.init_image_height)
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
    
    # def disable_curvefit_start(self):
    #     if self.dytran_use_t0_curvefit_checkbutton_var.get(): 
    #        self.dytran_curvefit_start_time_entry.config(state = tk.DISABLED)    
    #     else:
    #        self.dytran_curvefit_start_time_entry.config(state = tk.NORMAL)  
                
    
    def create_dytran_widgets(self):
        self.dytran_tab = tk.Frame(self.notebook)
        self.notebook.add(self.dytran_tab, text='Dytran')

        # dytran tab
        # First setup frames for dytran tab
        self.dytran_tab_left = tk.Frame(self.dytran_tab, bd=5, relief=tk.SUNKEN)
        self.dytran_tab_left.grid(row=0, column=0, sticky='NW')

        self.dytran_tab_right = tk.Frame(self.dytran_tab, bd=5, relief=tk.SUNKEN)
        self.dytran_tab_right.grid(row=0, column=1, sticky='NW')

        # Now within the dytran frames setup stuff
        self.load_dytran_checkbutton = tk.Checkbutton(
            self.dytran_tab_left,
            text='Load dytran',
            variable=self.load_dytran_checkbutton_var,
            command=self.load_fastDAQ_dytran)
        self.load_dytran_checkbutton.grid(row=0, column=0, columnspan=2, sticky='WE')

        self.dytran_label = tk.Label(self.dytran_tab_left, text='Dytran:')
        self.dytran_label.grid(row=1, column=0, sticky='WE')

        self.dytran_combobox = ttk.Combobox(self.dytran_tab_left, width=12)
        self.dytran_combobox.grid(row=1, column=1, sticky='WE')

        self.dytran_start_time_label = tk.Label(self.dytran_tab_left, text='Start Time:')
        self.dytran_start_time_label.grid(row=2, column=0, sticky='WE')

        self.dytran_start_time_entry = tk.Entry(self.dytran_tab_left, width=12)
        self.dytran_start_time_entry.insert(0, self.dytran_start_time)
        self.dytran_start_time_entry.grid(row=2, column=1, sticky='WE')

        self.dytran_ending_time_label = tk.Label(self.dytran_tab_left, text='Ending Time:')
        self.dytran_ending_time_label.grid(row=3, column=0, sticky='WE')

        self.dytran_ending_time_entry = tk.Entry(self.dytran_tab_left, width=12)
        self.dytran_ending_time_entry.insert(0, self.dytran_ending_time)
        self.dytran_ending_time_entry.grid(row=3, column=1, sticky='WE')
              
        self.dytran_timerange_checkbutton = tk.Checkbutton(
            self.dytran_tab_left, text='Full time window',
            variable=self.dytran_timerange_checkbutton_var,
            command=self.draw_fastDAQ_dytran)
        self.dytran_timerange_checkbutton.grid(row=4, column=0, columnspan=2, sticky='WE')

        self.dytran_lowpass_checkbutton = tk.Checkbutton(
            self.dytran_tab_left, text='Apply lowpass filter',
            variable=self.dytran_lowpass_checkbutton_var,
            command=self.draw_fastDAQ_dytran)
        self.dytran_lowpass_checkbutton.grid(row=5, column=0, columnspan=2, sticky='WE')

        self.dytran_plot_t0_checkbutton = tk.Checkbutton(
            self.dytran_tab_left,
            text='Show t0',
            variable=self.dytran_plot_t0_checkbutton_var,
            command=self.draw_fastDAQ_dytran)
        self.dytran_plot_t0_checkbutton.grid(row=6, column=0, columnspan=2, sticky='WE')
        
        self.dytran_plot_curvefit_checkbutton = tk.Checkbutton(
            self.dytran_tab_left,
            text='Show curve fit',
            variable=self.dytran_curvefit_checkbutton_var,
            command=self.draw_fastDAQ_dytran)
        self.dytran_plot_curvefit_checkbutton.grid(row=7, column=0, columnspan=2, sticky='WE')
        
        self.dytran_use_t0_curvefit_checkbutton = tk.Checkbutton(
            self.dytran_tab_left,
            text='Use t0 for fit start',
            state = tk.DISABLED,
            variable = self.dytran_use_t0_curvefit_checkbutton_var,
            command = self.draw_fastDAQ_dytran)
        self.dytran_use_t0_curvefit_checkbutton.grid(row=8, column=0, columnspan=2, sticky='WE')

        self.dytran_calc_t_end_curvefit_checkbutton = tk.Checkbutton(
            self.dytran_tab_left,
            text='Est. compression for fit end',
            state = tk.DISABLED,
            variable = self.dytran_calc_t_end_curvefit_checkbutton_var,
            command = self.draw_fastDAQ_dytran)
        self.dytran_calc_t_end_curvefit_checkbutton.grid(row=9, column=0, columnspan=2, sticky='WE')

        self.dytran_curvefit_start_time_label = tk.Label(self.dytran_tab_left, text='Curve fit t_start:')
        self.dytran_curvefit_start_time_label.grid(row= 10, column=0, sticky='WE')
        
        self.dytran_curvefit_start_time_entry = tk.Entry(self.dytran_tab_left, width=12)
        self.dytran_curvefit_start_time_entry.insert(0, self.dytran_curvefit_start_time)
        self.dytran_curvefit_start_time_entry.grid(row=10, column=1, sticky='WE')
        
        self.dytran_curvefit_window_time_label = tk.Label(self.dytran_tab_left, text='Fit window size:')
        self.dytran_curvefit_window_time_label.grid(row= 11, column=0, sticky='WE')
        
        self.dytran_curvefit_window_time_entry = tk.Entry(self.dytran_tab_left, width=12)
        self.dytran_curvefit_window_time_entry.insert(0, self.dytran_curvefit_window_time)
        self.dytran_curvefit_window_time_entry.grid(row=11, column=1, sticky='WE')

        self.reload_fastDAQ_dytran_button = tk.Button(self.dytran_tab_left, text='reload',
                                                      command=self.draw_fastDAQ_dytran)
        self.reload_fastDAQ_dytran_button.grid(row=12, column=0, sticky='WE')
