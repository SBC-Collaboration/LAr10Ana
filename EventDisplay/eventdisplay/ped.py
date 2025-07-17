#!/coupp/app/home/coupp/anaconda3/bin/python
'''pico-event-display

to run: python ped.py
may need to add to your paths:
 export PATH=/coupp/app/home/coupp/anaconda3/bin:$PATH
 export PYTHONPATH=/coupp/app/home/coupp/anaconda3/pkgs:$PYTHONPATH

'''

# Imports
import os
os.umask(6)
import re
import time
import getpass
import logging
import linecache
import matplotlib
import numpy as np
from pylab import *
from pathlib import Path
import tkinter as tk
from tkinter import ttk
import platform

matplotlib.use('TkAgg')
from tkinter import messagebox
from tkinter import filedialog
from PIL import PngImagePlugin
import tarfile
import zipfile
import sys
from tabs.camera import Camera
from tabs.piezo import Piezo
from tabs.logviewer import LogViewer
from tabs.configuration import Configuration
from tabs.analysis import Analysis
from tabs.three_d_bubble import ThreeDBubble
from tabs.scintillation import Scintillation
from GetEvent import GetEvent

try:
    from ctypes import windll
    windll.shcore.SetProcessDpiAwareness(1)
except:
    pass

PngImagePlugin.MAX_TEXT_CHUNK = 2000
# verbosity = logging.DEBUG
verbosity = logging.INFO
tar_postfixes = ['.tar', '.tar.gz', '.tgz']

# Defines message box for errors
class PopUpHandler(logging.Handler):
    def __init__(self):
        logging.Handler.__init__(self)
        self.setFormatter(logging.Formatter(fmt='%(message)s'))

    def emit(self, message):
        messagebox.showerror('error', self.format(message))


# Sets width/height/zoom settings for each tk frame (rectangular window where widgets can be placed)
class Application(Camera, Piezo, LogViewer, Configuration, Analysis, ThreeDBubble, Scintillation):
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)

        self.time_start = 0
        self.init_image_width = 400  # 659
        self.init_image_height = 625  # 494
        self.native_image_width = None
        self.native_image_height = None
        self.max_zoom = 3

        self.scbcanvas = tk.Canvas(ROOT, borderwidth=0)
        self.scbframe = tk.Frame(self.scbcanvas)
        self.scbframe.grid()
        self.vsb = tk.Scrollbar(ROOT, orient="vertical", command=self.scbcanvas.yview)
        self.hsb = tk.Scrollbar(ROOT, orient="horizontal", command=self.scbcanvas.xview)
        self.scbcanvas.configure(yscrollcommand=self.vsb.set)
        self.scbcanvas.configure(xscrollcommand=self.hsb.set)

        self.vsb.pack(side="right", fill="y")
        self.hsb.pack(side="bottom", fill="x")
        self.scbcanvas.pack(side="left", fill="both", expand=True)
        self.scbcanvas.create_window((4, 4), window=self.scbframe, anchor="nw", tags="self.scbframe")

        self.scbframe.bind("<Configure>", self.onFrameConfigure)
        self.scbframe.bind_all("<MouseWheel>", self._on_mousewheel)

        # Logger Initial Settings#
        self.logger = logging.getLogger('ped')
        self.formatter = logging.Formatter(fmt='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
                                           datefmt='%Y-%m-%d %H:%M:%S')
        self.console = logging.StreamHandler()
        self.console.setFormatter(self.formatter)
        self.logger.addHandler(self.console)
        self.logger.setLevel(verbosity)

        self.pop_up = PopUpHandler()
        self.pop_up.setLevel(logging.ERROR)
        self.logger.addHandler(self.pop_up)

        #################### USER-SPECIFIC-PATHS ###################
        # # for Data Directory if no config file found
        # # ex: self.raw_init_directory = 'C:\\Users\\User_name\\Research\\PICO\\PICO_DATA\\30l-16-data'

        self.raw_init_directory = ''
        self.zipped_event = None
        self.zip_flag = False

        # # Directory paths based on self.dataset
        # ped_directory: where the ped code is stored on the machine
        # scan_directory: where the handscans are stored
        # npy_directory: where the .npy files (raw_events.npy, reco_events.npy) and merged_all.txt
        # config_file_directory: where the configuration files are placed, should be a folder in the same directory as ped_directory
        # extraction_path: where tar files are unpacked
        self.dataset = 'SBC-25'
        self.reco_filename = 'reco_events.npy'
        self.ped_directory = os.path.split(os.path.split(os.path.abspath(__file__))[0])[0]
        self.scan_directory = os.path.join(self.ped_directory, 'scan_output_' + self.dataset)
        self.npy_directory = os.path.join(self.ped_directory, 'npy', self.dataset)
        self.config_file_directory = os.path.join(self.ped_directory, 'configs')
        self.extraction_path = os.path.join(self.ped_directory, 'scratch', getpass.getuser())
        self.init_run = ''

        # Errors will be appended to this string and displayed for each event
        self.error = ''

        # Default config values
        self.image_naming_conventions = ['cam0_image0.png', 'cam0image  0.bmp', 'cam1-img00.png']
        self.image_orientations = ['0', '90', '180', '270']
        self.plc_temp_var = 'TE32'
        self.images_relative_path = 'Images'
        self.image_naming_convention = self.image_naming_conventions[0]
        self.num_cams = 4
        self.image_orientation = self.image_orientations[3]
        self.first_frame = '30'
        self.init_frame = '50'
        self.last_frame = '70'
        self.piezo = 'Piezo7'
        self.ped_config_file_path_var = os.path.join(self.config_file_directory, self.dataset + '-ped_config.txt')

        # # try to load config file from self.dataset, then load user input self.raw_init_directory if no config file found
        self.load_config_values(self.ped_config_file_path_var)

        # Directory paths based on config file
        # raw_init_directory: where the data for an individual dataset is stored
        # base_directory: where the datasets are stored, created from raw_directory
        # log_directory: where the log files are placed, should be a folder in the raw_directory
        ##########################################################

        self.get_raw_events()

        self.source_button_var = tk.IntVar(value=-1)
        self.nbub_button_var = tk.IntVar(value=-1)
        self.do_handscan_checkbutton_var = tk.BooleanVar(value=False)
        self.use_cut_file_checkbutton_var = tk.BooleanVar(value=False)
        self.draw_crosshairs_var = tk.BooleanVar(value=False)
        self.invert_checkbutton_var = tk.BooleanVar(value=False)
        self.diff_checkbutton_var = tk.BooleanVar(value=False)
        self.antialias_checkbutton_var = tk.BooleanVar(value=False)
        self.load_dytran_checkbutton_var = tk.BooleanVar(value=False)
        self.piezo_plot_t0_checkbutton_var = tk.BooleanVar(value=False)
        self.load_fastDAQ_piezo_checkbutton_var = tk.BooleanVar(value=False)
        self.isgoodtrigger_checkbutton_var = tk.BooleanVar(value=True)
        self.crosshairsgood_checkbutton_var = tk.BooleanVar(value=True)
        self.load_initial_data_checkbutton_var = tk.BooleanVar(value=False)
        self.recons_signal_var = tk.BooleanVar(value=False)
        for i in range(9):
            self.grid_rowconfigure(i, weight=1)

        self.scanner_name = tk.StringVar()

        # PLC vars
        self.temp_label = tk.StringVar()

        # event_info.sbc vars
        self.trigger_type = -1
        self.trigger_type_label = tk.StringVar()
        self.pset_label = tk.StringVar()
        self.livetime_label = tk.StringVar()

        self.run = None
        self.event = None
        self.row_index = -1
        self.cuts = []
        self.selected_events = None
        self.selected_reco_indices = None
        self.reco_events = None
        self.reco_row = None

        # Initial Functions
        self.create_widgets()

        Camera.__init__(self)
        Piezo.__init__(self)
        LogViewer.__init__(self)
        Analysis.__init__(self)
        ThreeDBubble.__init__(self)
        Configuration.__init__(self)
        Scintillation.__init__(self)

        # self.load_reco()

        # Initial Functions
        self.initialize_widget_values()
        self.reset_event()

        #Icon setup for Windows, Mac and Linux
        try:
           iconName = 'PICO'
           windowSystem = self.master.tk.call("tk", "windowingsystem")
           if windowSystem == "x11" or windowSystem =='aqua': # Unix and Mac
               iconName = iconName + ".gif"
               iconImage = tk.PhotoImage(file=os.path.join(self.ped_directory, iconName))
               ROOT.call('wm', 'iconphoto', ROOT._w, iconImage)
           else: # Windows
               iconName += ".ico"
               ROOT.iconbitmap(os.path.join(self.ped_directory, iconName))
        except:
            print('Unable to add PICO icon to GUI')

    def onFrameConfigure(self, event):
        '''Reset the scroll region to encompass the inner frame'''
        self.scbcanvas.configure(scrollregion=self.scbcanvas.bbox("all"))

    def _on_mousewheel(self, event):
        self.scbcanvas.yview_scroll(int(-1*(event.delta/120)), "units")

    def click_tab(self, event):
        # Get Tab Clicked
        try:
            tab_clicked = self.notebook.tk.call(self.notebook._w, 'identify', 'tab', event.x, event.y)
        except tk.TclError:
            if event.widget.identify(event.x, event.y) == 'label':
                index = event.widget.index('@%d,%d' % (event.x, event.y))
                tab_clicked = event.widget.tab(index, 'text')
            else:
                tab_clicked = ''
        except Exception as e:
            print(e, 'Line: ', sys.exc_info()[2].tb_lineno)

        # Click on Log Viewer Tab
        if tab_clicked == 3 or tab_clicked == 'Log Viewer':
            # Show Log Viewer Widgets
            self.bottom_frame_5.grid(row=1, column=0, sticky='NW')
            self.fra2.pack(side='top')

            # Hide Other Bottom Frames
            self.bottom_frame_1.grid_remove()
            self.bottom_frame_2.grid_remove()
            self.bottom_frame_3.grid_remove()
            self.bottom_frame_4.grid_remove()
        elif tab_clicked == '':
            return
        else:
            # Show Other Bottom Frames
            self.bottom_frame_1.grid(row=1, column=0, sticky='NW')
            self.bottom_frame_2.grid(row=1, column=1, sticky='NW')
            self.bottom_frame_3.grid(row=1, column=2, sticky='NW')
            self.bottom_frame_4.grid(row=1, column=3, sticky='NW')

            # Hide Log Viewer Widgets
            self.bottom_frame_5.grid_remove()
            self.fra2.pack_forget()

    def initialize_widget_values(self):
        values = sorted(self.reco_events.dtype.names) if self.reco_events is not None else ('')
        self.add_display_var_combobox['values'] = values
        self.path_ped_config_directory.insert(0, self.config_file_directory)
        self.ped_config_file_path_combobox['values'] = self.get_configs()
        self.piezo_combobox['values'] = [self.piezo]
        self.piezo_combobox.current(0)
        # self.dytran_combobox['values'] = [self.dytran]
        # self.dytran_combobox.current(0)
        self.piezo_selector_combobox['values'] = [self.piezo]
        self.piezo_selector_combobox.current(0)
        if os.path.isfile(self.ped_config_file_path_var):
            self.ped_config_file_path_combobox.insert(0, os.path.basename(self.ped_config_file_path_var))

    def _resolve_path(self, raw):
        """
        Take relative path and resolve it to work with windows, macOS and Linux.
        """
        p = Path(os.path.expandvars(str(raw))).expanduser()

        if not p.is_absolute():
            p = Path(self.ped_directory) / p

        return p.resolve()

    # reads config file and sets given values, otherwise sets default values
    def load_config_values(self, path):
        values = [None] * 16

        defaults = []
        defaults.insert(0, str(self.raw_init_directory))
        defaults.insert(1, str(self.extraction_path))
        defaults.insert(2, str(self.init_run))
        defaults.insert(3, 'TE32')
        defaults.insert(4, 'Images')
        defaults.insert(5, 0)
        defaults.insert(6, 4)
        defaults.insert(7, 3)
        defaults.insert(8, '30')
        defaults.insert(9, '50')
        defaults.insert(10, '70')
        defaults.insert(11, 'Piezo0')
        defaults.insert(12, 'Dytran')
        defaults.insert(13, '500')
        defaults.insert(14, '1000')
        defaults.insert(15, '500')

        if os.path.isfile(path):
            f = open(path)
            counter = 0

            for line in f:
                parsed_line = line.rstrip('\n')
                if counter % 2 == 1:
                    values[int((counter - 1) / 2)] = parsed_line
                counter += 1
        else:
            self.logger.error(
                'Config file not found at: ' + self.ped_config_file_path_var + '\n\nDefault values set \n\nUSER-SPECIFIC-PATHS self.raw_init_directory set to: ' + self.raw_init_directory)

        for num in range(len(values)):
            if values[num] is None:
                values[num] = defaults[num]

        # Resolve base directories
        self.raw_init_directory = str(self._resolve_path(values[0]))
        self.raw_directory = self.raw_init_directory
        # Order matters here, we need to try removing -daqdata first, then -data
        self.dataset = os.path.basename(self.raw_init_directory).removesuffix('-daqdata').removesuffix('-data')

        # Relative directories
        self.scan_directory = str(self._resolve_path(f"scan_output_{self.dataset}"))
        self.npy_directory  = str(self._resolve_path(Path('npy') / self.dataset))
        self.log_directory  = str((Path(self.raw_init_directory) / 'logs').resolve())
        self.extraction_path = str((self._resolve_path(values[1]) / getpass.getuser()).resolve())

        # Scalar settings
        self.init_run = values[2]
        self.plc_temp_var = values[3]
        self.images_relative_path = values[4]
        self.image_naming_convention = self.image_naming_conventions[int(values[5])]
        self.num_cams = int(values[6])
        self.image_orientation = self.image_orientations[int(values[7])]
        self.first_frame = values[8]
        self.init_frame = values[9]
        self.last_frame = values[10]
        self.piezo = values[11]
        self.dytran = values[12]
        self.radius = values[13]
        self.positive_z = values[14]
        self.negative_z = values[15]
        self.frame = self.init_frame

        self.base_directory = os.path.dirname(self.raw_init_directory)


    def reset_event(self):
        self.run = ''
        self.reco_row = None
        self.row_index = -1
        self.zip_flag = False
        self.increment_event(1)

        # Change run to init_run, if starting run is set
        if self.init_run != '':
            self.load_run(self.init_run, 0)

    def get_raw_events(self):
        # print('raw path: ', self.raw_directory, '\n')
        user_date = '{}_{}'.format(getpass.getuser(), time.strftime('%a_%b_%d_%H_%M_%S_%Y'))
        try:
            self.raw_events = np.load(os.path.join(self.npy_directory, 'raw_events.npy'))
        except FileNotFoundError:
            try:
                os.system("python \"{}\" \"{}\" \"{}\"".format(os.path.join(self.ped_directory, "convert_raw_to_npy_run_by_run.py"), self.raw_directory, self.npy_directory))
                os.system("python \"{}\" \"{}\"".format(os.path.join(self.ped_directory, "merge_raw_run_npy.py"), self.npy_directory))
                # reco_path = os.path.join(self.npy_directory, self.reco_filename)
                # if not os.path.isfile(reco_path):
                #     reco_response = messagebox.askyesno('No raw_event.npy or reco_events.npy files found', 'NPY files being created now. \nWould you like to select a reco_events file?')
                #     if reco_response == 0:
                #         merged_all_response = messagebox.askyesno('Create a reco_events.npy file?', 'Would you like to select a merged_all file to create reco data?')
                #         if merged_all_response == 0:
                #             self.get_raw_events()
                #         else:
                #             merged_all = filedialog.askopenfilename(initialdir = self.npy_directory, title = "Select a merged_all File", filetypes = (("Text files", "*.txt*"), ("all files", "*.*")))
                #             os.system("python \"{}\" \"{}\" \"{}\" \"{}\" \"{}\"".format(os.path.join(self.ped_directory, "convert_reco_to_npy_and_reindex_raw_npy.py"), self.npy_directory, self.npy_directory, merged_all, user_date))
                #             self.reco_filename = 'reco_events_{}.npy'.format(user_date)
                #     else:
                #         reco_filename = filedialog.askopenfilename(initialdir = self.npy_directory, title = "Select a reco_events File", filetypes = (("NPY Files", "*.npy*"), ("all files", "*.*")))
                #         if os.path.isfile(reco_filename):
                #             reco_directory = os.path.split(reco_filename)[0]
                #             if os.path.normpath(reco_directory) == os.path.normpath(self.npy_directory):
                #                 self.reco_filename = os.path.split(reco_filename)[1]
                #             else:
                #                 directory_response = messagebox.askyesno('Selected File Not in NPY Directory', 'Current Dataset is: {}. \nSelected reco file is not in the NPY Directory for the current dataset. \nWould you still like to use the selected reco_events file?'.format(self.dataset))
                #                 if directory_response == 1:
                #                     self.reco_filename = os.path.split(reco_filename)[1]
                #                 else:
                #                     self.get_raw_events()
                #         else:
                #             self.get_raw_events()
                #         self.reco_filename = os.path.split(reco_filename)[1]
                # else:
                #     messagebox.showinfo('No raw_events.npy file found', 'NPY files being created now. \nreco_event.npy file found. Reco data will be loaded from npy directory')
            except FileNotFoundError:
                # this error should be handled when it crops up in the code
                raise FileNotFoundError
        try:
            self.raw_events = np.load(os.path.join(self.npy_directory, 'raw_events.npy'))
        except FileNotFoundError:
            # this error should be handled when it crops up in the code
            raise FileNotFoundError

    def convert_reco_from_merged_all(self):
        user_date = '{}_{}'.format(getpass.getuser(), time.strftime('%a_%b_%d_%H_%M_%S_%Y'))
        merged_all = filedialog.askopenfilename(
            initialdir=self.npy_directory,
            title="Select a merged_all File",
            filetypes=(("Text files", "*.txt*"), ("all files", "*.*")))
        if len(merged_all) > 0:
            os.system("python \"{}\" \"{}\" \"{}\" \"{}\" \"{}\"".format(
                os.path.join(self.ped_directory, "convert_reco_to_npy_and_reindex_raw_npy.py"),
                self.npy_directory,
                self.npy_directory,
                merged_all,
                user_date))
        else:
            return

        self.reco_filename = 'reco_events_{}.npy'.format(user_date)

        self.get_raw_events()
        self.update_num_cams()
        self.load_reco()
        self.reset_event()
        self.set_init_dataset()
        self.reset_images()
        self.reset_cuts()
        self.load_3d_bubble_data()

    def select_reco_file(self):
        if platform.system() == 'Darwin':       # macOS
            reco_filename = filedialog.askopenfilename(initialdir=self.npy_directory, title="Select a reco_events File")
        else:        
            reco_filename = filedialog.askopenfilename(initialdir=self.npy_directory, title="Select a reco_events File", filetypes=(("NPY Files", "*.npy*"), ("all files", "*.*")))

        if os.path.isfile(reco_filename):
            reco_directory = os.path.split(reco_filename)[0]
            if os.path.normpath(reco_directory) == os.path.normpath(self.npy_directory):
                self.reco_filename = os.path.split(reco_filename)[1]
            else:
                directory_response = messagebox.askyesno('Selected File Not in NPY Directory', 'Current Dataset is: {}. \nSelected reco file is not in the NPY Directory for the current dataset. \nWould you still like to use the selected reco_events file?'.format(self.dataset))
                if directory_response == 1:
                    self.reco_filename = os.path.split(reco_filename)[1]
                else:
                    return
        else:
            return
        # self.ped_directory = os.path.split(os.path.split(os.path.abspath(__file__))[0])[0]

        self.show_var.set(False)

        self.get_raw_events()
        self.load_reco()
        self.reset_event()
        self.set_init_dataset()
        self.reset_images()
        self.reset_cuts()
        self.load_3d_bubble_data()

    def load_event_sbc(self):
        selected = ["event_info"]
        self.path = os.path.join(self.raw_directory, self.run)

        try:
            event_info = GetEvent(self.path, self.event, *selected)["event_info"]
            livetime = event_info["ev_livetime"][0]
            pset = event_info["pset"][0]
            trigger_source = event_info["trigger_source"][0]

            self.trigger_type_label.set(f'trig: {trigger_source}')
            self.pset_label.set(f'pset: {pset:.1f}')
            self.livetime_label.set(f'lt: {livetime:.1f}')

        except:
            self.trigger_type_label.set('trigger: N/A')
            self.pset_label.set('pset: N/A')
            self.livetime_label.set('lt: N/A')
            self.error += 'cannot find event_info.sbc\n'

    def plc_text_zip_loader(self, path:str) -> None:
        with self.zipped_event.open(path) as file:
            try:
                fields = file.readline()
                fields = file.readline().split()
                fields = [field.decode() for field in fields]
                index = fields.index(self.plc_temp_var)
                entries = file.readline()
                entries = file.readline()
                entries = file.readline()
                entries = file.readline()
                entries = file.readline()
                entries = file.readline().split()
                entries = [entry.decode() for entry in entries]
                self.temp_label.set(self.plc_temp_var + ': {:.1f}'.format(float(entries[index])))
            except:
                self.error += 'cannot find ' + self.plc_temp_var + ' in PLC log file (via zip)\n'
                self.temp_label.set(self.plc_temp_var + ': N/A')

    def load_plc_text(self):
        return
        if self.zip_flag:
            path = os.path.join(self.run, str(self.event), 'PLClog.txt')
            self.plc_text_zip_loader(path)
        else:
            path = os.path.join(self.raw_directory, self.run, str(self.event), 'PLClog.txt')
            try:
                fields = linecache.getline(path, 2).split()
                index = fields.index(self.plc_temp_var)
                entries = linecache.getline(path, 7).split()
                self.temp_label.set(self.plc_temp_var + ': {:.1f}'.format(float(entries[index])))
                #print(f'fields: {fields}\n index: {index}\n entries:{entries}\n')
            except ValueError:
                self.temp_label.set(self.plc_temp_var + ': N/A')
                self.error += 'cannot find ' + self.plc_temp_var + ' in PLC log file\n'

    def archive_file_helper(self, run_path:str, extract_path:str) -> bool:
        """helper function to find and extract zip and tarfiles

        Args:
            run_path (str): where the runfile is
            extract_path (str): wher to extract archive contents usually
            the scratch directory

        Returns:
            bool: True if found archived file type
        """        
        test_path = run_path + '.zip'
        if os.path.exists(test_path):
            #self.logger.error('zip file found. Attempting to unzip...')
            self.zipped_event = zipfile.ZipFile(test_path, 'r')
            self.raw_directory = self.raw_directory
            self.zip_flag = True
            #self.logger.error('Zip flag has set to True')
            #zipfile.ZipFile(test_path).extractall(self.raw_directory )
            return True
        else:
            tar_postfixes = ['.tar', '.tar.gz', '.tgz']
            for ext in tar_postfixes:
                test_path = run_path + ext
                if os.path.exists(test_path):
                    self.logger.error('tar file found. Attempting to untar...')
                    t = tarfile.open(test_path, 'r')
                    self.raw_directory = self.extraction_path
                    t.extractall(self.raw_directory)
                    return True
        return False
                

    # Deal with possibility that run folders are tarred
    def handle_run_folder_format(self):
        self.zip_flag = False
        self.raw_directory = self.raw_init_directory
        run_folder_path = os.path.join(self.raw_directory, self.run)
        if not os.path.exists(run_folder_path):
            run_scratch_path = os.path.join(self.extraction_path, self.run, '0', 'Event.txt')
            if os.path.exists(run_scratch_path):
                self.raw_directory = self.extraction_path
                self.logger.info('Non-empty run folder found in scratch dir')
            else:
                #self.logger.error('Non-empty run folder does not exist. Searching for, and unzip/untarring  zip/tar version...')
                archive_file_found = False
                run_archive_path = os.path.join(self.raw_directory, self.run)
                self.logger.info(run_archive_path)
                #Check for tar/zip files
                archive_file_found = self.archive_file_helper(run_archive_path, self.extraction_path)
                if not archive_file_found:
                    self.logger.error('zip/tar file not found.')

    def reload_run(self):
        if self.selected_events is None:
            self.row_index = self.get_row(self.raw_events)
        else:
            try:
                self.row_index = self.get_row(self.selected_events)
            except IndexError:
                self.logger.error('disabling cuts: requested run does not satisfy them')
                self.selected_events = None
                self.selected_reco_indices = None
                self.row_index = self.get_row(self.raw_events)
        self.load_reco_row()

    def load_run(self, run, event):
        if run == self.run and event == self.event:
            self.logger.info('no action taken (run and event are unchanged)')
        else:
            if len(np.argwhere((self.raw_events['run'] == run) & (self.raw_events['ev'] == event))) == 0:
                self.logger.error('invalid request: run {}, event {} does not exist'.format(run, event))
                self.update_run_entry()
                return
            self.logger.info('going to run {}, event {}'.format(run, event))

            prevrun = self.run
            self.run = run
            if self.run != prevrun:
                self.handle_run_folder_format()
            self.event = event

            if self.selected_events is None:
                self.row_index = self.get_row(self.raw_events)
            else:
                try:
                    self.row_index = self.get_row(self.selected_events)
                except IndexError:
                    self.logger.error('disabling cuts: requested run does not satisfy them')
                    self.selected_events = None
                    self.selected_reco_indices = None
                    self.row_index = self.get_row(self.raw_events)
            self.update_run_entry()
            self.load_reco_row()
            if self.zip_flag:
                self.image_directory = os.path.join(run, str(event), self.images_relative_path)
            else:
                self.image_directory = os.path.join(self.raw_directory, run, str(event), self.images_relative_path)
            self.reset_images()

    def add_display_var(self, var):
        if (self.reco_events is not None) and (var not in self.reco_events.dtype.names):
            self.logger.error('requested variable not in reco data: ' + var)
            return

        if var in [label['text'] for label, text, value in self.display_vars]:
            return

        label = tk.Label(self.bottom_frame_2, text=var)
        label.grid(row=len(self.display_vars) + 2, column=0)

        text = tk.StringVar(value=self.reco_row[var]) if self.reco_row else tk.StringVar(value='N/A')

        value = tk.Label(self.bottom_frame_2, textvariable=text, width=8)
        value.grid(row=len(self.display_vars) + 2, column=1, sticky='W')

        self.display_vars.append((label, text, value))

    def add_cut(self):
        field = ttk.Combobox(self.bottom_frame_1, width=3, values=sorted(self.reco_events.dtype.names))
        field.insert(0, 'nbub')
        field.grid(row=8 + len(self.cuts), column=0, columnspan=2, sticky='WE')

        operator = ttk.Combobox(self.bottom_frame_1, width=3, values=('>', '>=', '==', '<=', '<', '!='))
        operator.insert(0, '>=')
        operator.grid(row=8 + len(self.cuts), column=2, sticky='WE')

        value = tk.Entry(self.bottom_frame_1, width=5)
        value.insert(0, '0')
        value.grid(row=8 + len(self.cuts), column=3, sticky='WE')

        self.cuts.append((field, operator, value))

    def remove_cut(self):
        if not self.cuts:
            return

        for widget in self.cuts.pop():
            widget.destroy()

        self.apply_cuts()
        if self.selected_events == None:
            self.reload_run()

    def remove_all_cuts(self):
        if not self.cuts:
            return

        while self.cuts:
            for widget in self.cuts.pop():
                widget.destroy()

        self.apply_cuts()
        if self.selected_events == None:
            self.reload_run()

    def reset_cuts(self):
        for field, operator, value in self.cuts:
            field.delete(0, tk.END)
            operator.delete(0, tk.END)
            value.delete(0, tk.END)

        self.apply_cuts()
        if self.selected_events == None:
            self.reload_run()

    def apply_cuts(self):
        self.selected_events = None
        self.selected_reco_indices = None

        if self.reco_events is None:
            self.logger.error('cannot apply cuts, reco data not found')
            return

        selection = []
        for field, operator, value in self.cuts:
            if field.get() == '' and operator.get() == '' and value.get() == '':
                continue
            if field.get() not in self.reco_events.dtype.names:
                self.logger.error('requested variable not in reco data')
                field.delete(0, tk.END)
                return
            dtype = self.reco_events[field.get()].dtype.str
            selection.append('(self.reco_events["{}"] {} {})'.format(
                field.get(),
                operator.get(),
                repr(value.get()) if 'U' in dtype else value.get()))  # add quotes if field datatype is string

        if len(selection) > 0:
            selection = eval(' & '.join(selection))
            selected_event_indices = np.where(selection)
            if len(selected_event_indices) > 0:
                selected_event_indices = selected_event_indices[0]
                self.selected_events = self.reco_events[selected_event_indices]
            else:
                self.logger.error('no events pass cuts')
                self.reset_cuts()
                return

            # exec('self.selected_events = self.reco_events[{}]'.format(' & '.join(selection)))
            _, unique_rows = np.unique(self.selected_events[['run', 'ev']], return_index=True)
            self.selected_events = self.selected_events[unique_rows]  # get rid of multiple nbub entries
            self.selected_reco_indices = selected_event_indices[unique_rows]
            # print('unique rows:', unique_rows)
            # print('reco indices for unique rows:', self.selected_reco_indices)

            row = self.get_row(self.raw_events)
            prevrun = self.run
            try:
                events_left = self.raw_events[['run', 'ev']][row:]
                run, event = np.intersect1d(self.selected_events[['run', 'ev']], events_left)[0]
            except IndexError:
                self.logger.error('reached final event: starting over1')
                run, event = self.selected_events[['run', 'ev']][0]
            self.run = run
            if self.run != prevrun:
                self.handle_run_folder_format()
            self.event = event
            self.reco_row = None
            self.row_index = self.get_row(self.selected_events) - 1
            self.increment_event(1)
        else:
            self.selected_events = None
            self.selected_reco_indices = None
            # self.row_index = self.get_row(self.raw_events)
            self.row_index = 0

    def add_file_cut(self):
        self.cut_file_label = tk.Label(self.bottom_frame_1, text='Select .txt file from npy directory')
        self.cut_file_label.grid(row=8 + len(self.cuts), column=0, columnspan=2, sticky='WE')

        self.cut_file_combobox = ttk.Combobox(self.bottom_frame_1, width=3)
        self.cut_file_combobox['values'] = self.get_cut_files()
        self.cut_file_combobox.insert(0, 'No File Selected')
        self.cut_file_combobox.grid(row=8 + len(self.cuts), column=2, columnspan=3, sticky='WE')

        self.cuts.append((self.cut_file_label, self.cut_file_combobox))

    # Returns a list of all config files in the config directory
    def get_cut_files(self):
        all_cut_files = os.listdir(self.npy_directory)
        cut_files = []
        for file in all_cut_files:
            fileRegex = re.compile(r'.*\.txt$')
            if fileRegex.match(file):
                cut_files.append(str(file))
        return cut_files

    def apply_file_cuts(self):
        self.selected_events = None
        self.selected_reco_indices = None

        if self.reco_events is None:
            self.logger.error('cannot apply cuts, reco data not found')
            return

        selection = []

        self.cut_file = os.path.join(self.npy_directory, self.cut_file_combobox.get())

        try:
            cuts = np.loadtxt(self.cut_file, dtype=str, ndmin=2)
        except ValueError:
            self.logger.error('Too many fields, no cuts applied')
            return

        for cut in cuts:
            selection.append('((self.reco_events["{}"] {} "{}") & (self.reco_events["{}"] {} {}))'.format('run', '==', cut[0], 'ev', '==', cut[1]))

        if len(selection) > 0:
            selection = eval(' | '.join(selection))
            selected_event_indices = np.where(selection)
            if len(selected_event_indices) > 0:
                selected_event_indices = selected_event_indices[0]
                self.selected_events = self.reco_events[selected_event_indices]
            else:
                self.logger.error('no events pass cuts')
                self.reset_cuts()
                return

            # exec('self.selected_events = self.reco_events[{}]'.format(' & '.join(selection)))
            _, unique_rows = np.unique(self.selected_events[['run', 'ev']], return_index=True)
            self.selected_events = self.selected_events[unique_rows]  # get rid of multiple nbub entries
            self.selected_reco_indices = selected_event_indices[unique_rows]
            # print('unique rows:', unique_rows)
            # print('reco indices for unique rows:', self.selected_reco_indices)

            row = self.get_row(self.raw_events)
            prevrun = self.run
            try:
                events_left = self.raw_events[['run', 'ev']][row:]
                run, event = np.intersect1d(self.selected_events[['run', 'ev']], events_left)[0]
            except IndexError:
                self.logger.error('reached final event: starting over2')
                run, event = self.selected_events[['run', 'ev']][0]
            self.run = run
            if self.run != prevrun:
                self.handle_run_folder_format()
            self.event = event
            self.reco_row = None
            self.row_index = self.get_row(self.selected_events) - 1
            self.increment_event(1)
        else:
            self.selected_events = None
            self.selected_reco_indices = None
            # self.row_index = self.get_row(self.raw_events)
            self.row_index = 0

    def get_row(self, events):
        rows = np.argwhere((events['run'] == self.run) & (events['ev'] == self.event)).ravel()
        if len(rows) > 0:
            return rows[0]
        else:
            return -1

    # For moving forward and backwards through events, loads the next appropriate row of data
    def load_reco_row(self, ibub=None):
        # print('In Load reco row')
        # print('self.run: ' + self.run)
        # print('self.ev: ' + str(self.event))
        # print('row_index: ', str(self.row_index))
        if self.reco_events is not None:
            self.reco_availability_label.config(text='Loaded')
        else:
            self.reco_availability_label.config(text='Not Loaded')
            return

        if self.selected_events is not None:
            # here the row index is the index relative to the selected event list
            self.reco_row = self.selected_events[self.row_index]
        else:
            # here the row index is the index relative to the raw event list
            reco_index_fast = self.raw_events[self.row_index]['reco index']
            reco_index = self.get_row(self.reco_events)    # should give same result as above, but slower
            if reco_index_fast != reco_index:
                print('fast reco index is wrong: ')
                print(' slow reco index: ', reco_index)
                print(' fast reco index: ', reco_index_fast)
            if reco_index >= 0:
                self.reco_row = self.reco_events[reco_index]
            else:
                self.reco_row = None
                self.toggle_reco_widgets(state=tk.DISABLED)
                for _, text, _ in self.display_vars:
                    text.set('N/A')
                return
            # print('reco_index: ', str(reco_index))

        if ibub:
            # print('  reco.run: ', self.reco_row['run'])
            # print('  reco.ev: ', self.reco_row['ev'])
            # print('  nbub: ', self.reco_row['nbub'])
            offset = ibub - 1 if ibub > 1 else 0
            if self.selected_reco_indices is not None:
                # print('sel reco indices: ', self.selected_reco_indices)
                row_fast_selection = self.selected_reco_indices[self.row_index]
            else:
                row_fast_selection = self.raw_events[self.row_index]['reco index']
            row = self.get_row(self.reco_events)     # should give same result as above, but slower
            if row_fast_selection != row:
                print('fast row index is wrong: ')
                print(' slow row index: ', row)
                print(' fast row selection index: ', row_fast_selection)
                print('   raw row index: ', self.row_index)
                print('   ibub: ', ibub)
                print('   offset: ', offset)
            self.reco_row = self.reco_events[row + int(offset)]
            # print('  ->reco.run: ', self.reco_row['run'])
            # print('  ->reco.ev: ', self.reco_row['ev'])
            # print('  ->nbub: ', self.reco_row['nbub'])

        self.toggle_reco_widgets(state=tk.NORMAL)
        for label, text, _ in self.display_vars:
            var = label['text']
            try:
                dtype = self.reco_row[var].dtype.str
                text.set('{:.4f}'.format(self.reco_row[var]) if 'f' in dtype else self.reco_row[var])
            except:
                text.set('N/A')

    def use_cut_file(self):
        if self.use_cut_file_checkbutton_var.get():
            self.add_cut_button['state'] = tk.DISABLED
            self.apply_cuts_button['state'] = tk.DISABLED
            self.reset_cuts_button['state'] = tk.DISABLED
            self.remove_cut_button['state'] = tk.DISABLED
            self.add_file_cut_button['state'] = tk.NORMAL
            self.apply_file_cuts_button['state'] = tk.NORMAL
        else:
            self.add_cut_button['state'] = tk.NORMAL
            self.apply_cuts_button['state'] = tk.NORMAL
            self.reset_cuts_button['state'] = tk.NORMAL
            self.remove_cut_button['state'] = tk.NORMAL
            self.add_file_cut_button['state'] = tk.DISABLED
            self.apply_file_cuts_button['state'] = tk.DISABLED

        self.remove_all_cuts()

    def increment_event(self, step):
        self.error = ''
        if self.selected_events is None:
            events = self.raw_events
        else:
            events = self.selected_events

        if (self.row_index + step) < 0:
            self.logger.error('reached first event: stopping here')
            self.reset_event()
            return

        if (self.row_index + step) >= len(events):
            self.logger.error('reached final event: starting over3')
            self.reset_event()
            return

        self.row_index += step

        prevrun = self.run
        self.run = events[self.row_index]['run']
        if self.run != prevrun:
            self.handle_run_folder_format()
        self.event = events[self.row_index]['ev']

        self.update_run_entry()
        if self.zip_flag:
            self.image_directory = os.path.join(self.run, str(self.event), self.images_relative_path)
        else:
            self.image_directory = os.path.join(self.raw_directory, self.run, str(self.event), self.images_relative_path)

        self.load_reco_row()
        self.reset_images()

        self.load_3d_bubble_data()

        # printing out all error messages for an event at once
        if not self.error == '':
            self.logger.error('This event had the following errors:\n' + self.error)

    ##find all dataset folders in the base directory of the form X*(size)-XX(year)-data
    def get_datasets(self):
        files = os.listdir(self.base_directory)
        i = 0
        while i < files.__len__():
            fileRegex = re.compile('\\w*-\\w*-data')
            if not fileRegex.match(files[i]):
                files.remove(files[i])
                i = i - 1
            i = i + 1
        return files

    # Returns a list of all config files in the config directory
    def get_configs(self):
        all_files = os.listdir(self.config_file_directory)
        files = []
        for file in all_files:
            fileRegex = re.compile(r'.*config(?!.*_EXAMPLE\.txt).*\.txt$')
            if fileRegex.match(file):
                files.append(str(file))
        return files

    # Returns a list of all custom reco files in the npy directory
    def get_custom_reco(self):
        all_files = os.listdir(self.npy_directory)
        files = []
        for file in all_files:
            fileRegex = re.compile(r'.*reco_events.*\.npy$')
            if fileRegex.match(file):
                files.append(str(file))
        return files

    def update_run_entry(self):
        self.run_entry.delete(0, tk.END)
        self.run_entry.insert(0, self.run)
        self.event_entry.delete(0, tk.END)
        self.event_entry.insert(0, self.event)

    def set_init_dataset(self):
        counter = 0
        for dataset in self.get_datasets():
            if self.raw_directory.endswith(dataset):
                self.dataset_select.current(counter)
            counter += 1

    def load_reco(self):
        self.reco_row = None
        self.reco_events = None

        path = os.path.join(self.npy_directory, self.reco_filename)
        if not os.path.isfile(path):
            self.logger.error('cannot find {}, reco data will be disabled'.format(self.reco_filename))
            self.reco_availability_label.config(text='Not Loaded')
            self.toggle_reco_widgets(state=tk.DISABLED)
            for _, text, _ in self.display_vars:
                text.set('N/A')
            return

        self.logger.info('using reco data from {}'.format(path))

        events = np.load(path)
        if len(events) == 0:
            self.logger.error('could not find raw data for any reco events')
            return

        self.reco_events = events

    def do_handscan(self):
        if not os.path.exists(self.scan_directory):
            os.mkdir(self.scan_directory)

        state = tk.NORMAL
        if self.do_handscan_checkbutton_var.get():
            file = 'scan_{}_{}_{}.txt'.format(self.run, getpass.getuser(), time.strftime('%a_%b_%d_%H_%M_%S_%Y'))
            self.scan_file = os.path.join(self.scan_directory, file)
            self.scanner_name.set('scanner: ' + getpass.getuser())
            messagebox.showinfo(message='Beginning handscanning. Hitting submit will automatically save to a text file move you to the next event. Make sure you have no active cuts on reco data unless that is intended!')
        else:
            if os.path.isfile(self.scan_file) and (os.stat(self.scan_file).st_size == 0):
                os.remove(self.scan_file)
            state = tk.DISABLED
            self.scanner_name.set('')

        for i in range(0, 8):
            self.nbub_button[i].config(state=state)

        for button in self.source_buttons:
            button.config(state=state)

        self.nbub_button_var.set(-1)
        self.source_button_var.set(-1)
        self.isgoodtrigger_checkbutton_var.set(True)
        self.crosshairsgood_checkbutton_var.set(True)
        self.comment_entry.delete(0, tk.END)
        self.comment_entry.insert(0, '')

        self.isgoodtrigger_button.config(state=state)
        self.crosshairsgood_button.config(state=state)
        self.comment_label.config(state=state)
        self.comment_entry.config(state=state)
        self.submit_scan_button.config(state=state)

    def submit_scan(self):
        if ((str(self.scanner_name.get())[9:] == '') or
                (str(self.source_button_var.get()) == '-1') or
                (str(self.nbub_button_var.get()) == '-1')):
            self.logger.error('did not complete scan selections')
            return

        with(open(self.scan_file, 'a+')) as file:
            file.seek(0)
            if not file.read(1):
                file.write('Output from ped hand scanning.\n')
                file.write(
                    'run  ev  scanner  scan_source  scan_nbub  scan_trigger  scan_crosshairsgood  scan_comment\n')
                file.write('%s  %d  %s  %d  %d  %d  %d  %s\n1\n\n\n')

            file.write(self.run + '  ' +
                       str(self.event) + '  ' +
                       str(self.scanner_name.get())[9:] + '  ' +
                       str(self.source_button_var.get()) + ' ' +
                       str(self.nbub_button_var.get()) + ' ' +
                       str(int(self.isgoodtrigger_checkbutton_var.get())) + ' ' +
                       str(int(self.crosshairsgood_checkbutton_var.get())) + ' ' +
                       '\'' + self.comment_entry.get() + '\'\n')

        self.nbub_button_var.set(-1)
        self.source_button_var.set(-1)
        self.isgoodtrigger_checkbutton_var.set(True)
        self.crosshairsgood_checkbutton_var.set(True)
        self.comment_entry.delete(0, tk.END)
        self.comment_entry.insert(0, '')

        self.increment_event(1)

    def create_widgets(self):
        self.notebook = ttk.Notebook(self.scbframe, padding=[0, 0, 0, 0])

        # Function to Run When Clicking Tabs
        self.notebook.bind('<Button-1>', self.click_tab)

        self.notebook.grid(row=0, column=0, columnspan=5)

        # Setup frames to be used on the bottom

        self.bottom_frame_1 = tk.Frame(self.scbframe, bd=5, relief=tk.SUNKEN)
        self.bottom_frame_1.grid(row=1, column=0, sticky='NW')
        self.bottom_frame_2 = tk.Frame(self.scbframe, bd=5, relief=tk.SUNKEN)
        self.bottom_frame_2.grid(row=1, column=1, sticky='NW')
        self.bottom_frame_3 = tk.Frame(self.scbframe, bd=5, relief=tk.SUNKEN)
        self.bottom_frame_3.grid(row=1, column=2, sticky='NW')
        self.bottom_frame_4 = tk.Frame(self.scbframe, bd=5, relief=tk.SUNKEN)
        self.bottom_frame_4.grid(row=1, column=3, sticky='NW')
        self.bottom_frame_5 = tk.Frame(self.scbframe, bd=5, relief=tk.SUNKEN)  # Frame for Log Navigation
        self.bottom_frame_5.grid(row=1, column=0, sticky='NW')  # Frame for Log Navigation

        self.run_label = tk.Label(self.bottom_frame_1, text='run:')
        self.run_label.grid(row=1, column=0, sticky='WE')

        self.run_entry = tk.Entry(self.bottom_frame_1, width=12)
        self.run_entry.grid(row=1, column=1, sticky='WE')

        self.event_label = tk.Label(self.bottom_frame_1, text='event:')
        self.event_label.grid(row=1, column=2, sticky='WE')

        self.event_entry = tk.Entry(self.bottom_frame_1, width=5)
        self.event_entry.grid(row=1, column=3, sticky='WE')

        self.go_button = tk.Button(self.bottom_frame_1, text='Go', command=self.load_run)
        self.go_button['command'] = lambda: self.load_run(self.run_entry.get(), int(self.event_entry.get()))
        self.go_button.grid(row=1, column=4, sticky='WE')

        self.back_event = tk.Button(self.bottom_frame_1, text='back event', command=lambda: self.increment_event(-1))
        self.back_event.grid(row=2, column=0, columnspan=2, sticky='WE')

        self.forward_event = tk.Button(self.bottom_frame_1, text='forward event',
                                       command=lambda: self.increment_event(1))
        self.forward_event.grid(row=2, column=2, columnspan=2, sticky='WE')

        self.back_1000events_button = tk.Button(self.bottom_frame_1, text='back 1000 events')
        self.back_1000events_button['command'] = lambda: self.increment_event(-1000)
        self.back_1000events_button.grid(row=3, column=0, columnspan=2, sticky='WE')

        self.forward_1000events_button = tk.Button(self.bottom_frame_1, text='forward 1000 events')
        self.forward_1000events_button['command'] = lambda: self.increment_event(1000)
        self.forward_1000events_button.grid(row=3, column=2, columnspan=2, sticky='WE')

        self.fill_trigger_type = tk.Label(self.bottom_frame_1, textvariable=self.trigger_type_label, width=11)
        self.fill_trigger_type.grid(row=4, column=0, sticky='WE')

        self.fill_pset = tk.Label(self.bottom_frame_1, textvariable=self.pset_label, width=10)
        self.fill_pset.grid(row=4, column=1, sticky='WE')

        self.fill_temp = tk.Label(self.bottom_frame_1, textvariable=self.temp_label, width=10)
        self.fill_temp.grid(row=4, column=2, sticky='WE')

        self.fill_livetime = tk.Label(self.bottom_frame_1, textvariable=self.livetime_label, width=10)
        self.fill_livetime.grid(row=4, column=3, sticky='WE')

        self.reset_cuts_button = tk.Button(self.bottom_frame_1, text='reset cuts', command=self.reset_cuts)
        self.reset_cuts_button.grid(row=6, column=0, columnspan=2, sticky='WE')

        self.apply_cuts_button = tk.Button(self.bottom_frame_1, text='apply cuts', command=self.apply_cuts)
        self.apply_cuts_button.grid(row=6, column=2, columnspan=2, sticky='WE')

        self.add_cut_button = tk.Button(self.bottom_frame_1, text='add cut', command=self.add_cut)
        self.add_cut_button.grid(row=5, column=0, columnspan=2, sticky='WE')

        self.remove_cut_button = tk.Button(self.bottom_frame_1, text='delete cut', command=self.remove_cut)
        self.remove_cut_button.grid(row=5, column=2, columnspan=2, sticky='WE')

        self.add_file_cut_button = tk.Button(self.bottom_frame_1, text='add cut from .txt', width=25, command=self.add_file_cut)
        self.add_file_cut_button['state'] = tk.DISABLED
        self.add_file_cut_button.grid(row=7, column=0, columnspan=2, sticky='WE')

        self.apply_file_cuts_button = tk.Button(self.bottom_frame_1, text='apply cuts from .txt', command=self.apply_file_cuts)
        self.apply_file_cuts_button['state'] = tk.DISABLED
        self.apply_file_cuts_button.grid(row=7, column=2, columnspan=2, sticky='WE')

        self.use_cut_file_checkbutton = tk.Checkbutton(self.bottom_frame_1,
            text='Use Cut File?',
            variable=self.use_cut_file_checkbutton_var,
            command=self.use_cut_file)
        self.use_cut_file_checkbutton.grid(row=7, column=4, sticky='WE')

        self.display_reco_label = tk.Label(self.bottom_frame_2, text='Variables from merged_all')
        self.display_reco_label.grid(row=0, column=0, sticky='WE')

        self.add_display_var_combobox = ttk.Combobox(self.bottom_frame_2)
        self.add_display_var_combobox.grid(row=1, column=0, sticky='WE')

        self.add_display_var_button = tk.Button(
            self.bottom_frame_2,
            text='add',
            command=lambda: self.add_display_var(self.add_display_var_combobox.get()))
        self.add_display_var_button.grid(row=1, column=1, sticky='WE')

        self.display_vars = []
        self.add_display_var('nbub')
        #     self.add_display_var('getBub_success')
        self.add_display_var('fastDAQ_t0')
        self.add_display_var('te')

        self.reco_label = tk.Label(self.bottom_frame_2, text='Reco status:')
        self.reco_label.grid(row=0, column=2, sticky='WE')

        self.reco_availability_label = tk.Label(self.bottom_frame_2, text='Not Loaded')
        self.reco_availability_label.grid(row=1, column=2, sticky='WE')

        self.select_reco_button = tk.Button(self.bottom_frame_2, text='Select reco file', command=self.select_reco_file)
        self.select_reco_button.grid(row=2, column=2, sticky='WE')

        self.convert_reco_button = tk.Button(self.bottom_frame_2, text='Create reco file', command=self.convert_reco_from_merged_all)
        self.convert_reco_button.grid(row=3, column=2, sticky='WE')

        self.back_frame_button = tk.Button(self.bottom_frame_3, text='back frame')
        self.back_frame_button['command'] = lambda: self.load_frame(int(self.frame) - 1)
        self.back_frame_button.grid(row=0, column=0, sticky='WE')

        self.forward_frame_button = tk.Button(self.bottom_frame_3, text='forward frame')
        self.forward_frame_button['command'] = lambda: self.load_frame(int(self.frame) + 1)
        self.forward_frame_button.grid(row=0, column=1, sticky='WE')

        self.reset_images_button = tk.Button(self.bottom_frame_3, text='reset image', command=self.reset_images)
        self.reset_images_button.grid(row=0, column=2, sticky='WE')

        self.first_frame_button = tk.Button(self.bottom_frame_3, text='first frame')
        self.first_frame_button['command'] = lambda: self.load_frame(self.first_frame)
        self.first_frame_button.grid(row=1, column=0, sticky='WE')

        self.last_frame_button = tk.Button(self.bottom_frame_3, text='last frame')
        self.last_frame_button['command'] = lambda: self.load_frame(self.last_frame)
        self.last_frame_button.grid(row=1, column=1, sticky='WE')

        self.trig_frame_button = tk.Button(self.bottom_frame_3, text='trig frame')
        self.trig_frame_button['command'] = lambda: self.load_frame(self.init_frame)
        self.trig_frame_button.grid(row=1, column=2, sticky='WE')

        self.antialias_checkbutton = tk.Checkbutton(self.bottom_frame_3,
            text='antialias',
            variable=self.antialias_checkbutton_var,
            command=self.update_images)
        self.antialias_checkbutton.grid(row=2, column=0, sticky='WE')

        self.diff_checkbutton = tk.Checkbutton(self.bottom_frame_3,
            text='diff frame',
            variable=self.diff_checkbutton_var,
            command=self.update_images)
        self.diff_checkbutton.grid(row=2, column=1, sticky='WE')

        self.invert_checkbutton = tk.Checkbutton(self.bottom_frame_3,
            text='invert',
            variable=self.invert_checkbutton_var,
            command=self.update_images)
        self.invert_checkbutton.grid(row=2, column=2, sticky='WE')

        self.draw_crosshairs_button = tk.Checkbutton(self.bottom_frame_3,
            text='draw crosshairs',
            variable=self.draw_crosshairs_var,
            command=self.draw_crosshairs,
            state=tk.DISABLED)
        self.draw_crosshairs_button.grid(row=3, column=0, sticky='WE')

        self.make_video_button = tk.Button(self.bottom_frame_3,
            text='make video',
            command=self.make_video)
        self.make_video_button.grid(row=4, column=0, sticky='WE')

        self.make_video_label = tk.Label(self.bottom_frame_3, text="cam")
        self.make_video_label.grid(row=4, column=1, sticky='WE')

        self.make_video_entry = tk.Entry(self.bottom_frame_3, width=10)
        self.make_video_entry.insert(0, '0')
        self.make_video_entry.grid(row=4, column=2, sticky='WE')

        self.do_handscan_checkbutton = tk.Checkbutton(self.bottom_frame_4,
            text='do handscan',
            variable=self.do_handscan_checkbutton_var,
            command=self.do_handscan)
        self.do_handscan_checkbutton.grid(row=0, column=0, columnspan=2, sticky='WE')

        self.scanner_name_label = tk.Label(self.bottom_frame_4, textvariable=self.scanner_name)
        self.scanner_name_label.grid(row=0, column=2, columnspan=2, sticky='WE')

        self.nbub_button = []
        max_bub = 8
        for i, text in zip(range(0, max_bub), [' bubbles'] * (max_bub - 1) + ['+ bubbles']):
            self.nbub_button.append(
                tk.Radiobutton(
                    self.bottom_frame_4,
                    text=str(i) + text,
                    state=tk.DISABLED,
                    variable=self.nbub_button_var,
                    command=self.change_nbub,
                    value=i))
            self.nbub_button[i].grid(row=1 + i, column=0, columnspan=2, sticky='WE')

        self.source_buttons = []
        for i, text in enumerate(['bulk event', 'wall', 'dome', 'bellows region', 'other']):
            button = tk.Radiobutton(self.bottom_frame_4, state=tk.DISABLED, variable=self.source_button_var, value=i)
            button['text'] = text
            button.grid(row=i + 1, column=2, columnspan=2, sticky='WE')
            self.source_buttons.append(button)

        self.isgoodtrigger_button = tk.Checkbutton(self.bottom_frame_4, variable=self.isgoodtrigger_checkbutton_var)
        self.isgoodtrigger_button['text'] = 'Is good trigger?'
        self.isgoodtrigger_button['state'] = tk.DISABLED
        self.isgoodtrigger_button.grid(row=7, column=2, columnspan=2, sticky='WE')

        self.crosshairsgood_button = tk.Checkbutton(self.bottom_frame_4, variable=self.crosshairsgood_checkbutton_var)
        self.crosshairsgood_button['text'] = 'Crosshairs good?'
        self.crosshairsgood_button['state'] = tk.DISABLED
        self.crosshairsgood_button.grid(row=8, column=2, columnspan=2, sticky='WE')

        self.comment_label = tk.Label(self.bottom_frame_4, text='Comment:', state=tk.DISABLED)
        self.comment_label.grid(row=9, column=0, sticky='WE')

        self.comment_entry = tk.Entry(self.bottom_frame_4, width=15, state=tk.DISABLED)
        self.comment_entry.insert(0, '')
        self.comment_entry.grid(row=9, column=1, columnspan=2, sticky='WE')

        self.submit_scan_button = tk.Button(self.bottom_frame_4, state=tk.DISABLED, text='Submit and Go')
        self.submit_scan_button['command'] = self.submit_scan
        self.submit_scan_button.grid(row=9, column=3, sticky='WE')

def on_closing():
    plt.close()
    ROOT.destroy()  # Close Window
    sys.exit()  # Stop Running Script


# Create Window and Set On Top
ROOT = tk.Tk()
ROOT.lift()
ROOT.attributes('-topmost', True)
ROOT.after_idle(ROOT.attributes, '-topmost', False)

# Size and Position of Window
WIDTH = 1650
HEIGHT = 975
X = 0
Y = 0

# Update Window and Run
ROOT.geometry('%dx%d+%d+%d' % (WIDTH, HEIGHT, X, Y))
ROOT.title('PED')
# ROOT.iconbitmap(os.path.join('..', 'PICO.ico'))
APP = Application(ROOT)
ROOT.protocol("WM_DELETE_WINDOW", on_closing)
APP.mainloop()
