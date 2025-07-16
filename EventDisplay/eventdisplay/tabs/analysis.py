import tkinter as tk
from pandas import DataFrame
import gc
import os
import numpy as np
import time
import matplotlib.pyplot as plt
from time import sleep
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from tkinter import ttk, messagebox, DISABLED, NORMAL
from matplotlib.figure import Figure
from tkinter import *
import sys
from PIL import Image, ImageTk
import scipy.signal
import matplotlib
import matplotlib.pyplot as plt
# matplotlib.use('TkAgg')
matplotlib.use('Agg')
from scipy.signal import butter, sosfilt
import getpass
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
from GetEvent import GetEvent

class Analysis(tk.Frame):
     def __init__(self, master=None):
          tk.Frame.__init__(self, master)


          self.last_hline = 0
          self.freq_cutoff_low = 10000
          self.freq_cutoff_high = 100000

          self.piezo_analysis_beginning_time = 0
          self.piezo_analysis_ending_time = 1

          self.timerange_checkbutton_var = tk.BooleanVar(value=False)
          self.plot_t0_checkbutton_var = tk.BooleanVar(value=False)
          self.denoise_selector_var = tk.BooleanVar(value=False)
          self.increment_piezo_event = False

          self.create_analysis_widgets()
          self.analysis_canvas_setup()

          self.t1_start = 0
          self.flag1 = False
         
     def load_fastDAQ_analysis(self):
          if not self.load_initial_data_checkbutton_var.get():
               self.analysis_tab_graph.grid_forget()
               return
          else:
               self.analysis_tab_graph.grid(row=0, column=1, sticky='NW')

          if self.zip_flag:
               path = os.path.join(self.raw_directory, self.run, '.zip')

          path = os.path.join(self.raw_directory, self.run)

          try:
               selected = ["run_control", "acoustics"]
               self.fastDAQ_event = GetEvent(path, self.event, *selected)

               self.piezo_selector_combobox['values'] = [f"Channel {i+1}" for i in range(self.fastDAQ_event['acoustics']['Waveform'].shape[1])]

               self.check_t0_exist()
               self.draw_fastDAQ_analysis()
          except:
              print("not found")

          # Garbage Collecting
          gc.collect()

     def draw_fastDAQ_analysis(self):
          self.load_DAQ()
          self.scale_selector_setup()
          self.fill_t0_entry()
          self.reload_command()
          
     def piezo_selection(self):
         if not self.load_initial_data_checkbutton_var.get():
           messagebox.showerror('Error','Please load initial data.') #If initial data not loaded then a messagebox emerges 
           return
         
         self.piezo = self.piezo_selector_combobox.get()
         
         if self.piezo not in self.piezo_selector_combobox['values']:
           if self.piezo == '':
             messagebox.showerror('Error','Please make a piezo selection.')  # If the combobox in empty a messegebox emerges 
             self.piezo_selector_combobox.focus()
             return
           else:
             messagebox.showerror('Error','{} does not exist in this dataset.\nPlease change selection.'.format(self.piezo))  # If the piezo does not exit in the dataset a messegebox emerges 
             self.piezo_selector_combobox.focus()           
             return
             
         self.load_DAQ(self.piezo)

     def load_DAQ(self, piezo = None):
         if not self.load_initial_data_checkbutton_var.get(): 
           return        
         if self.piezo_selector_combobox.current() == -1:
             self.piezo_selector_combobox.current(0)
         self.f = self.fastDAQ_event['acoustics']['Waveform'][0][self.piezo_selector_combobox.current()]
         self.t = np.arange(len(self.f)) * (1 / self.fastDAQ_event['acoustics']['sample_rate'])
         self.piezo_selector_combobox['values'] = [f"Channel {i+1}" for i in range(self.fastDAQ_event['acoustics']['Waveform'].shape[1])]
         # self.f = 0.5*np.sin(2*np.pi*60*self.t) + np.sin(2*np.pi*1000*self.t)*(self.t>0) + 2.0 * np.sin(2*np.pi*2000*self.t)*(self.t>0) # for testing
         self.fft_variables()
         
     def fft_variables(self):
          self.n = len(self.t)                                                                 # How many data point there are
          self.dt = (self.t[-1] - self.t[0]) / ( self.n - 1 )
          self.fhat = np.fft.fft(self.f,self.n)                                                # Compute the FFT
          self.PSD = ( self.fhat * np.conj(self.fhat) / self.n ).real                          # Power spectrum density (power per frequency)
          self.freq = (1/(self.dt*self.n)) * np.arange(self.n)                                 # Creating x-axis of frequencies
          self.L = np.arange(1, np.floor(self.n/2), np.floor(self.n/2)/4000, dtype = 'int')    # Only plot the first half of the signal 
          self.PSD_max = np.max(self.PSD)
          self.PSD_min = np.min(self.PSD)

     def analysis_canvas_setup(self):
          # Create Figures, Axes, and Canvases
          plt.rcParams.update({'font.size': 8})
          plt.rcParams.update({'figure.autolayout': True})
          plt.rc('axes', titlesize=10)
          plt.rc('axes', labelsize=10)

          user = getpass.getuser()
          self.tmp_dir = os.path.join(self.ped_directory, "scratch", user, "tmp")

          # Create full directory path if it doesn't exist
          os.makedirs(self.tmp_dir, exist_ok=True)
                    
          image = Image.open(os.path.join(self.ped_directory, 'notfound.jpeg'))
          if image.mode == 'RGBA': # Image in RGBA, Pillow cannot deal with that, so convert to RGB
               r, g, b, a = image.split()
               image = Image.merge('RGB', (r, g, b))

          self.canvas11 = tk.Canvas(self.analysis_tab_graph, width=400, height=300)
          self.canvas11.image = self.canvas11.create_image(0, 0, anchor=tk.NW, image=None)
          self.canvas11.grid(row=0, column=1, sticky='NW')
          self.canvas11.photo = ImageTk.PhotoImage(image)
          self.canvas11.itemconfig(self.canvas11.image, image=self.canvas11.photo)

          self.canvas12 = tk.Canvas(self.analysis_tab_graph, width=400, height=300)
          self.canvas12.image = self.canvas12.create_image(0, 0, anchor=tk.NW, image=None)
          self.canvas12.grid(row=0, column=2, sticky='NW')
          self.canvas12.photo = ImageTk.PhotoImage(image)
          self.canvas12.itemconfig(self.canvas12.image, image=self.canvas12.photo)

          self.canvas21 = tk.Canvas(self.analysis_tab_graph, width=400, height=300)
          self.canvas21.image = self.canvas21.create_image(0, 0, anchor=tk.NW, image=None)
          self.canvas21.grid(row=1, column=1, sticky='NW')
          self.canvas21.photo = ImageTk.PhotoImage(image)
          self.canvas21.itemconfig(self.canvas21.image, image=self.canvas21.photo)

          self.canvas22 = tk.Canvas(self.analysis_tab_graph, width=400, height=300)
          self.canvas22.image = self.canvas22.create_image(0, 0, anchor=tk.NW, image=None)
          self.canvas22.grid(row=1, column=2, sticky='NW')
          self.canvas22.photo = ImageTk.PhotoImage(image)
          self.canvas22.itemconfig(self.canvas22.image, image=self.canvas22.photo)
          
          self.canvas31 = tk.Canvas(self.analysis_tab_graph, width=400, height=300)
          self.canvas31.image = self.canvas31.create_image(0, 0, anchor=tk.NW, image=None)
          self.canvas31.grid(row=0, column=3, sticky='NW')
          self.canvas31.photo = ImageTk.PhotoImage(image)
          self.canvas31.itemconfig(self.canvas31.image, image=self.canvas31.photo)
          
          self.canvas32 = tk.Canvas(self.analysis_tab_graph, width=400, height=300)
          self.canvas32.image = self.canvas32.create_image(0, 0, anchor=tk.NW, image=None)
          self.canvas32.grid(row=1, column=3, sticky='NW')
          self.canvas32.photo = ImageTk.PhotoImage(image)
          self.canvas32.itemconfig(self.canvas32.image, image=self.canvas32.photo)
        
     def plot_raw_and_filtered_signal(self):
          self.draw_raw_signal()
          self.plot_filtered_signal()
          
     def draw_raw_signal(self):
          if not self.load_initial_data_checkbutton_var.get():
               self.analysis_tab_graph.grid_forget()
               return
          else:          
               self.analysis_tab_graph.grid(row=0, column=1, sticky='NW')

               fig = plt.figure(figsize=(4, 3), dpi=100)
               ax = fig.add_subplot(111)

               ax.set_title('Raw Signal (Time Window)' + " " + str(self.run) + " " + str(self.event))
               ax.set_xlabel('[s]')
               ax.set_ylabel('[V]')
        

               self.piezo_analysis_beginning_time = float(self.piezo_analysis_beginning_time_entry.get())
               self.piezo_analysis_ending_time = float(self.piezo_analysis_ending_time_entry.get())

               self.piezo_analysis_ending_time_entry['state'] = tk.NORMAL
               self.piezo_analysis_beginning_time_entry['state'] = tk.NORMAL
               window = (self.t > self.piezo_analysis_beginning_time) & (self.t < self.piezo_analysis_ending_time)
               t = self.t[window]
               f = self.f[window]
               ax.plot(t[::int(len(t)/400)], f[::int(len(t)/400)])
               ax.margins(x=0)
                         
               if self.plot_t0_checkbutton_var.get():

                    try:             
                         if self.reco_row:
                              ax.axvline(x=self.reco_row['fastDAQ_t0'], linestyle='dashed', color='r', label='t0')
                              self.increment_piezo_event = True
                    except  ValueError:
                         return
                    try:
                         t0 = self.reco_row['fastDAQ_t0']             
                         if self.reco_row:
                              ax.axvline(x=t0, linestyle='dashed', color='r', label='t0')
                              self.increment_piezo_event = True
                         else:

                              if self.increment_piezo_event:
                                   self.error += 't0 unavailable: no reco data found for current event.'
                              else:
                                   self.logger.error('t0 unavailable: no reco data found for current event.')
                              self.plot_t0_checkbutton_var.set(False)
                              self.increment_piezo_event = False
                    except ValueError:
                         return


               fig.savefig(os.path.join(self.tmp_dir, 'raw.png'), dpi=100, bbox_inches='tight')
               plt.close(fig)
               
               image = Image.open(os.path.join(self.tmp_dir, 'raw.png'))
               if image.mode == 'RGBA': # Image in RGBA, Pillow cannot deal with that, so convert to RGB
                    r, g, b, a = image.split()
                    image = Image.merge('RGB', (r, g, b))

               self.canvas11.photo = ImageTk.PhotoImage(image)
               self.canvas11.itemconfig(self.canvas11.image, image=self.canvas11.photo)   
               
#                self.split_PSD()
               
     def log_convert(self, x):
          return 10**(x/10)
        
     def inv_log_convert(self, x): 
          return 10*np.log10(x)

         # Using FFT                    
     def plot_PSD(self, parameter = None):  
          if not self.load_initial_data_checkbutton_var.get():   
               self.analysis_tab_graph.grid_forget()
               return

          fig = plt.figure(figsize=(4, 3), dpi=100)
          ax = fig.add_subplot(111)
          
          ax.set_title('PSD')
          ax.set_xlabel('[KHz]')
          # ax.set_ylabel('[Power]')
               
          ax.plot(0.001 * self.freq[self.L],self.PSD[self.L])
          ax.margins(x=0)
          # print(len(self.freq[self.L]))
          ax.set_xscale('log')
          ax.set_xlim(left=1)
          ax.set_yscale('log')
          # ax.axhline(y= self.log_convert(self.slider.get()), color='r', linestyle='--')
          ax.axvline(x= 0.001 * float(self.freq_cutoff_low_entry.get()), linestyle='--', color='r', label='t0')
          ax.axvline(x= 0.001 * float(self.freq_cutoff_high_entry.get()), linestyle='--', color='r', label='t0')
                    
          fig.savefig(os.path.join(self.tmp_dir, 'psd.png'), dpi=100, bbox_inches='tight')
          plt.close(fig)
          
          image = Image.open(os.path.join(self.tmp_dir, 'psd.png'))
          if image.mode == 'RGBA': # Image in RGBA, Pillow cannot deal with that, so convert to RGB
               r, g, b, a = image.split()
               image = Image.merge('RGB', (r, g, b))

          self.canvas12.photo = ImageTk.PhotoImage(image)
          self.canvas12.itemconfig(self.canvas12.image, image=self.canvas12.photo)

     def denoise_signal_butter(self):
        # Use of PSD to filter out noise
        indices = self.PSD > self.log_convert(self.slider.get())      # Find all freqs with large power
        self.PSDclean = self.PSD * indices                            # Zero out all other
        fhat = indices * self.fhat                                         # Zero out small Fourier coeffs. in Y
        self.ffilt = ( np.fft.ifft(fhat) ).real                       # Inserse FFT for filtered time signal  
     
        # Use butter to filter frequency
        lowcut = float(self.freq_cutoff_low_entry.get())
        highcut = float(self.freq_cutoff_high_entry.get())
        nyq = len(self.ffilt)/(self.t[-1]-self.t[0])/2.0
        low = lowcut / nyq
        high = highcut / nyq
        if low <= 0:
            low = 0.0001
        if high >= 1:
            high = 0.9999
        sos = butter(6, [low, high], analog=False, btype='band', output='sos')
        self.ffilt = sosfilt(sos, self.ffilt)        

     def denoise_signal_fft(self):
        # Use of PSD to filter out freq
        indices_low_high = (self.freq > int(self.freq_cutoff_low_entry.get())) & (self.freq < int(self.freq_cutoff_high_entry.get()))   # Find frequencies outside the limits
        self.PSDclean = self.PSD * indices_low_high                                                                                     # Zero out freq outside limits
        fhat = indices_low_high * self.fhat                                                                                      
                                                                                   
        # Use of PSD to filter out noise
        indices = self.PSD > self.log_convert(self.slider.get())      # Find all freqs with large power
        self.PSDclean = self.PSDclean * indices                       # Zero out all other
        fhat = indices * fhat                                         # Zero out small Fourier coeffs. in Y
        self.ffilt = ( np.fft.ifft(fhat) ).real                       # Inserse FFT for filtered time signal

     def plot_filtered_signal(self):
          if not self.load_initial_data_checkbutton_var.get():
               self.analysis_tab_graph.grid_forget()
               return
          
          self.selected_filter = self.denoise_selector.get()     
          if self.selected_filter == 'Butter':
            self.denoise_signal_butter()
          else:
            self.denoise_signal_fft()
        
          self.piezo_analysis_beginning_time = float(self.piezo_analysis_beginning_time_entry.get())
          self.piezo_analysis_ending_time = float(self.piezo_analysis_ending_time_entry.get())

          fig = plt.figure(figsize=(4, 3), dpi=100)
          ax = fig.add_subplot(111)
        
          ax.set_title('Filtered Signal (Time Window)')
          ax.set_xlabel('[s]')
          ax.set_ylabel('[V]')
        

          self.piezo_analysis_ending_time_entry['state'] = tk.NORMAL
          self.piezo_analysis_beginning_time_entry['state'] = tk.NORMAL
          window = (self.t > self.piezo_analysis_beginning_time) & (self.t < self.piezo_analysis_ending_time)
          t = self.t[window]
          ffilt = self.ffilt[window]

               
          n = len(t)
          ax.plot(t[::int(n/40)], ffilt[::int(n/40)])
          ax.margins(x=0)

          if self.plot_t0_checkbutton_var.get():             
               try:
                    if self.reco_row:
                         ax.axvline(x=self.reco_row['fastDAQ_t0'], linestyle='dashed', color='r', label='t0')
                         self.increment_piezo_event = True
               except ValueError:
                    if self.increment_piezo_event:
                         self.error += 't0 unavailable: no reco data found for current event.'
                    else:
                         self.logger.error('t0 unavailable: no reco data found for current event.')
                    self.plot_t0_checkbutton_var.set(False)
                    self.increment_piezo_event = False

          fig.savefig(os.path.join(self.tmp_dir, 'filtered.png'), dpi=100, bbox_inches='tight')
          plt.close(fig)
          
          image = Image.open(os.path.join(self.tmp_dir, 'filtered.png'))
          if image.mode == 'RGBA': # Image in RGBA, Pillow cannot deal with that, so convert to RGB
               r, g, b, a = image.split()
               image = Image.merge('RGB', (r, g, b))
               
          self.canvas21.photo = ImageTk.PhotoImage(image)
          self.canvas21.itemconfig(self.canvas21.image, image=self.canvas21.photo)
              
     def plot_PSD_vs_time(self):
          if not self.load_initial_data_checkbutton_var.get(): 
               self.analysis_tab_graph.grid_forget()
               return

          PSD_list = []
          time_list = []
          for i in range(0, len(self.t), int(0.005 * len(self.t))):
               window = int(0.01 * len(self.t))
               f_ = self.ffilt[i:i+window]
               time_ = self.t[i:i+window]
               fhat_ = np.fft.fft(f_, len(time_))
               PSD_ = ( fhat_ * np.conj(fhat_) / len(time_) ).real            
               average_PSD = sum(PSD_)/window
               average_time = 0.5 * (time_[0] + time_[-1])
               PSD_list.append(average_PSD)
               time_list.append(average_time)

          fig = plt.figure(figsize=(4, 3), dpi=100)
          ax = fig.add_subplot(111)

          ax.clear()
          ax.set_title('Filtered PSD vs time')
          ax.set_xlabel('[s]')
          # ax.set_ylabel('[PSD]')
          ax.set_yscale('log')
          
          ax.plot(time_list, PSD_list)
          ax.margins(x=0)
          
          fig.savefig(os.path.join(self.tmp_dir, 'psd_vs_time.png'), dpi=100, bbox_inches='tight')
          plt.close(fig)
                         
          image = Image.open(os.path.join(self.tmp_dir, 'psd_vs_time.png'))
          if image.mode == 'RGBA': # Image in RGBA, Pillow cannot deal with that, so convert to RGB
               r, g, b, a = image.split()
               image = Image.merge('RGB', (r, g, b))
               
          self.canvas22.photo = ImageTk.PhotoImage(image)
          self.canvas22.itemconfig(self.canvas22.image, image=self.canvas22.photo)
          
     def split_PSD(self):
          if not self.load_initial_data_checkbutton_var.get():
            self.analysis_tab_graph.grid_forget()
            return
          else:          
            self.analysis_tab_graph.grid(row=0, column=1, sticky='NW')

          if self.post_t0_manual_input.get() == '':
               time_t0 = 0.00
          else:
               time_t0 = float(self.post_t0_manual_input.get())

          split_PSD_time_window = float(self.split_PSD_time_window_input.get())
          time_cut_pre = self.t <= self.t[0] + split_PSD_time_window
          time_cut_post = (self.t > time_t0) & (self.t < time_t0 + split_PSD_time_window)
          self.popup_PSD(self.f[time_cut_pre], self.t[time_cut_pre] , self.f[time_cut_post], self.t[time_cut_post])
          self.split_on_t0_signal(self.f[time_cut_pre], self.t[time_cut_pre] , self.f[time_cut_post], self.t[time_cut_post])
          
     def fill_t0_entry(self):
          try:
               t0 = self.reco_row['fastDAQ_t0']
               if self.reco_row:
                    self.post_t0_manual_input.delete(0, 'end')
                    self.post_t0_manual_input.insert(0, t0)

                    self.plot_t0_checkbutton.config(state=NORMAL)
          except (ValueError, TypeError):
             self.plot_t0_checkbutton.config(state=DISABLED)

          except (ValueError, TypeError):

             return
          
             
     def split_on_t0_signal(self, f1, t1, f2, t2):
        if not self.load_initial_data_checkbutton_var.get(): 
           self.analysis_tab_graph.grid_forget()
           return
           
        fig = plt.figure(figsize=(4, 3), dpi=100)
        ax = fig.add_subplot(111)


        ax.set_title('Raw Signal (Full)')
        ax.set_xlabel('[s]')
        ax.set_ylabel('[V]')    
        
        ax.plot(self.t, self.f)
        ax.margins(x=0)
                 
        f=f1
        t=t1
        ax.plot(t[::int(self.n/400)], f[::int(self.n/400)], color = 'purple', label = 'pre t0' ) 
        
        f=f2
        t=t2
        ax.plot(t[::int(self.n/400)], f[::int(self.n/400)], color = 'green', label = 'post t0' )  
        
         
        fig.savefig(os.path.join(self.tmp_dir, 'split_signal.png'), dpi=100, bbox_inches='tight')
        plt.close(fig)
          
        image = Image.open(os.path.join(self.tmp_dir, 'split_signal.png'))
        if image.mode == 'RGBA': # Image in RGBA, Pillow cannot deal with that, so convert to RGB
               r, g, b, a = image.split()
               image = Image.merge('RGB', (r, g, b))
               
        self.canvas31.photo = ImageTk.PhotoImage(image)
        self.canvas31.itemconfig(self.canvas31.image, image=self.canvas31.photo)
 
     
     def popup_PSD(self, f1, t1, f2, t2):          
          fig = plt.figure(figsize=(4, 3), dpi=100)
          ax = fig.add_subplot(111)
   
          ax.set_title('PSD comparison (pre/post t0)')
          ax.set_xlabel('[KHz]')
          #ax.set_ylabel('[V]')         

          f = f1
          t = t1
          n = len(t)                                                                 # How many data point there are
          dt = (t[-1] - t[0]) / ( n - 1 )
          fhat = np.fft.fft(f,n)                                                # Compute the FFT
          PSD = ( fhat * np.conj(fhat) / n ).real                          # Power spectrum density (power per frequency)
          freq = (1/(dt*n)) * np.arange(n)                                 # Creating x-axis of frequencies
          L = np.arange(1, np.floor(n/2), np.floor(n/2)/4000, dtype = 'int')    # Only plot the first half of the signal 
          PSD_max = np.max(PSD)
          PSD_min = np.min(PSD)          
                   
          ax.plot(0.001 * freq[L], scipy.signal.savgol_filter(PSD[L], 11, 1), color = 'purple', label = 'pre t0')

          f = f2
          t = t2
          n = len(t)                                                                 # How many data point there are
          dt = (t[-1] - t[0]) / ( n - 1 )
          fhat = np.fft.fft(f,n)                                                # Compute the FFT
          PSD = ( fhat * np.conj(fhat) / n ).real                          # Power spectrum density (power per frequency)
          freq = (1/(dt*n)) * np.arange(n)                                 # Creating x-axis of frequencies
          L = np.arange(1, np.floor(n/2), np.floor(n/2)/4000, dtype = 'int')    # Only plot the first half of the signal 
          PSD_max = np.max(PSD)
          PSD_min = np.min(PSD)
          
          ax.plot(0.001 * freq[L], scipy.signal.savgol_filter(PSD[L], 11, 1), color = 'green', label = 'post t0')
          ax.margins(x=0)
          ax.legend()
          
          ax.set_xscale('log')
          ax.set_xlim(left=1)
          ax.set_yscale('log') 
          
          fig.savefig(os.path.join(self.tmp_dir, 'split_psd.png'), dpi=100, bbox_inches='tight')
          plt.close(fig)
          
          image = Image.open(os.path.join(self.tmp_dir, 'split_psd.png'))
          if image.mode == 'RGBA': # Image in RGBA, Pillow cannot deal with that, so convert to RGB
               r, g, b, a = image.split()
               image = Image.merge('RGB', (r, g, b))
               
          self.canvas32.photo = ImageTk.PhotoImage(image)
          self.canvas32.itemconfig(self.canvas32.image, image=self.canvas32.photo)
          
          
     def Sc_notation_label(self):
       self.label.config(text=f"{10**(self.slider.get()/10):.2e}")
         
     def update_PSDvsTime_FiltedSignal(self, parameter = None):
          self.plot_PSD()
          self.plot_filtered_signal()
          self.plot_PSD_vs_time()
          self.Sc_notation_label()
          
     def update(self, parameter = None):
          y = self.slider.get()
          ymax = self.inv_log_convert(float(self.PSD_max))
          ymin = self.inv_log_convert(float(self.PSD_min))
          y_scaled = 30 + 210*(1-(y-ymin)/(ymax-ymin))
          self.canvas12.delete( self.last_hline )
          self.last_hline = self.canvas12.create_line(50, y_scaled, 380, y_scaled, dash=(9,3), fill='red', tag='hline')

          # self.plot_PSD()
          # if time.time() - self.t1_start < 1.0 or self.flag1 == False: 
          #      self.flag1 = True
          #      return
          # else:
          #      self.t1_start = time.time()
          #      self.plot_filtered_signal()
          #      self.plot_PSD_vs_time()
          #      self.Sc_notation_label()
                   
     def scale_selector_setup(self):
      if not self.load_initial_data_checkbutton_var.get():   #Wouldn't work correctly without these lines
            self.scale_tab.grid_forget()
            self.flag1 = False                    
            return
      else:      
        self.scale_tab.grid(row = 0, column = 3, sticky = 'NW', padx = 10)
        self.power_label = tk.Label(self.scale_tab, text='Power Cutoff')
        self.power_label.grid(row = 0, column = 0)
        self.label = tk.Label(self.scale_tab, text="1.00")
        self.label.grid(row = 1, column = 0)
        self.slider = Scale(self.scale_tab, from_= self.inv_log_convert(float(self.PSD_max)), to= self.inv_log_convert(float(self.PSD_min)), length= 550, command = self.update, activebackground = "grey65", troughcolor = 'white', showvalue = 0)
        self.slider.grid(row = 2, column = 0)
        self.slider.set((self.inv_log_convert(float(self.PSD_max)) + self.inv_log_convert(float(self.PSD_min)))/2)
        self.slider.bind("<ButtonRelease-1>", self.update_PSDvsTime_FiltedSignal)
        self.Sc_notation_label()
     
     def check_t0_exist(self):
          try:
               t0 = self.reco_row['fastDAQ_t0']
               self.plot_t0_checkbutton.config(state=NORMAL)
               return True
          except:
               self.plot_t0_checkbutton.config(state=DISABLED)
               return False

                                
     def reload_command(self):
          self.piezo_selection()
          if self.piezo not in self.piezo_selector_combobox['values']:    # If the combobox is empty then is not useful to run again all the plotting
               return
          self.check_t0_exist()
          self.draw_raw_signal()
          self.plot_PSD()
          self.plot_filtered_signal()
          self.plot_PSD_vs_time()
          self.split_PSD()
               
     def create_analysis_widgets(self):
        self.analysis_tab = tk.Frame(self.notebook)
        self.notebook.add(self.analysis_tab, text='Analysis')
         
        # setup frames within tab
        self.analysis_tab_main = tk.Frame(self.analysis_tab, bd=5, relief=tk.SUNKEN)
        self.analysis_tab_main.grid(row = 0, column = 0, sticky='NW')
        
        self.analysis_tab_graph = tk.Frame(self.analysis_tab, bd=5, relief=tk.SUNKEN)
        self.analysis_tab_graph.grid(row = 0, column = 1, sticky='NW')
        
        self.scale_tab = tk.Frame(self.analysis_tab, bd=5, relief=tk.SUNKEN)
        self.scale_tab.grid(row = 0, column = 2, sticky='NW')
       
        # Buttons & Labels        
        self.load_noiseData_checkbutton = tk.Checkbutton( self.analysis_tab_main, text='Load fastdaq', variable=self.load_initial_data_checkbutton_var, command=self.load_fastDAQ_analysis)
        self.load_noiseData_checkbutton.grid(row=0, column=0, columnspan=2, sticky='WE', padx=10)
        
        self.piezo_label = tk.Label(self.analysis_tab_main, text='Piezo:')
        self.piezo_label.grid(row = 1, column = 0, columnspan=2, sticky='WE')
        
        self.piezo_selector_combobox = ttk.Combobox(self.analysis_tab_main, width=12)
        self.piezo_selector_combobox.grid(row = 2, column = 0, columnspan=2, padx = 10)
        
        ttk.Separator(self.analysis_tab_main).grid(row = 3, column = 0, columnspan=2, sticky="ew", pady = 5)
        
        self.freq_cutoff_label = tk.Label(self.analysis_tab_main, text='Freq cutoff:')
        self.freq_cutoff_label.grid(row = 4, column = 0, columnspan=2, sticky = 'NEWS', padx = 10)

        self.freq_cutoff_low_label = tk.Label(self.analysis_tab_main, text='Low:', width=3)
        self.freq_cutoff_low_label.grid(row = 5, column = 0, sticky = 'E', padx = 10)

        self.freq_cutoff_low_entry = tk.Entry(self.analysis_tab_main, width=8)
        self.freq_cutoff_low_entry.insert(0, self.freq_cutoff_low)
        self.freq_cutoff_low_entry.grid(row = 5, column = 1, sticky='NEWS', padx = 10)

        self.freq_cutoff_high_label = tk.Label(self.analysis_tab_main, text='High:', width=3)
        self.freq_cutoff_high_label.grid(row = 6, column = 0, sticky = 'E', padx = 10)

        self.freq_cutoff_high_entry = tk.Entry(self.analysis_tab_main, width=8)
        self.freq_cutoff_high_entry.insert(0, self.freq_cutoff_high)
        self.freq_cutoff_high_entry.grid(row = 6, column = 1, sticky='NEWS', padx = 10)
        
        self.denoise_selector_label = tk.Label(self.analysis_tab_main, text='Filter type:')
        self.denoise_selector_label.grid(row= 7, column=0, columnspan=2, sticky='NEWS', padx = 10)

        self.denoise_selector = ttk.Combobox(self.analysis_tab_main, values= ['FFT', 'Butter'] , width=12)
        self.denoise_selector.grid(row=8, column=0, columnspan=2, sticky='NEWS', padx = 10)
        self.denoise_selector.current(1)
         
        ttk.Separator(self.analysis_tab_main).grid(row = 9, column = 0, columnspan=2, sticky="ew", pady = 5)
        
        self.piezo_analysis_time_label = tk.Label(self.analysis_tab_main, text='View window:')
        self.piezo_analysis_time_label.grid(row=10, column=0, columnspan=2, sticky='NEWS')

        self.piezo_analysis_beginning_time_label = tk.Label(self.analysis_tab_main, text='Begin:', width=4)
        self.piezo_analysis_beginning_time_label.grid(row=11, column=0, sticky='E')

        self.piezo_analysis_beginning_time_entry = tk.Entry(self.analysis_tab_main, width=7)
        self.piezo_analysis_beginning_time_entry.insert(0, self.piezo_analysis_beginning_time)
        self.piezo_analysis_beginning_time_entry.grid(row=11, column=1, sticky='NEWS', padx = 10)

        self.piezo_analysis_ending_time_label = tk.Label(self.analysis_tab_main, text='End:', width=4)
        self.piezo_analysis_ending_time_label.grid(row=12, column=0, sticky='E')

        self.piezo_analysis_ending_time_entry = tk.Entry(self.analysis_tab_main, width=4)
        self.piezo_analysis_ending_time_entry.insert(0, self.piezo_analysis_ending_time)
        self.piezo_analysis_ending_time_entry.grid(row=12, column=1, sticky='NEWS', padx = 10)

        ttk.Separator(self.analysis_tab_main).grid(row = 13, column = 0, columnspan=2, sticky="ew", pady = 5)
                
        self.split_PSD_time_window_input_label = tk.Label(self.analysis_tab_main, text='PSD comparison:')
        self.split_PSD_time_window_input_label.grid(row = 14, column = 0, columnspan=2, sticky = 'NEWS', padx = 10)

        self.split_PSD_time_window_input_label2 = tk.Label(self.analysis_tab_main, text='t window:', width=7)
        self.split_PSD_time_window_input_label2.grid(row = 15, column = 0, sticky = 'E', padx = (10,0))

        self.split_PSD_time_window_input = tk.Entry(self.analysis_tab_main, width=3) 
        self.split_PSD_time_window_input.insert(0, '0.01')
        self.split_PSD_time_window_input.grid(row = 15, column = 1, sticky = 'NEWS', padx = 10)
          
        self.post_t0_manual_input_label = tk.Label(self.analysis_tab_main, text='t0:', width=2) 
        self.post_t0_manual_input_label.grid(row = 16, column = 0, sticky = 'E', padx = 10)

        self.post_t0_manual_input = tk.Entry(self.analysis_tab_main, width=9)
        self.post_t0_manual_input.insert(0, '0.02')
        self.post_t0_manual_input.grid(row = 16, column = 1, sticky = 'NEWS', padx = 10)
        
        ttk.Separator(self.analysis_tab_main).grid(row = 17, column = 0, columnspan=2, sticky="ew", pady = 5)
        
        self.plot_t0_checkbutton = tk.Checkbutton(self.analysis_tab_main, text='Show t0', variable=self.plot_t0_checkbutton_var, command=self.plot_raw_and_filtered_signal)
        self.plot_t0_checkbutton.grid(row=18, column=0, columnspan=2, sticky='WE', padx=10)
        
        ttk.Separator(self.analysis_tab_main).grid(row = 19, column = 0, columnspan=2, sticky="ew", pady = 5)

        self.reload_button = tk.Button(self.analysis_tab_main,text='reload', command=self.reload_command)
        self.reload_button.grid(row = 20, column = 0, columnspan=2, sticky='NEWS', padx = 15, pady = 5)
