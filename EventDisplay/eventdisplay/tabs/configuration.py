# Imports
import os, getpass, time
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from io import open
class Configuration(tk.Frame):
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)

        # For config window
        self.coupp_checkbutton_var = tk.BooleanVar(value=False)
        self.use_custom_reco_checkbutton_var = tk.BooleanVar(value=False)

        self.raw_init_directory_path_var = tk.StringVar()
        self.scan_directory_path_var = tk.StringVar()
        self.npy_directory_path_var = tk.StringVar()

        self.raw_init_directory_path_var.set(self.raw_init_directory)
        self.scan_directory_path_var.set(self.scan_directory)
        self.npy_directory_path_var.set(self.npy_directory)

        # Initial Functions
        self.create_configuration_widgets()
        
        self.dataset_combobox_update()

    # Initial value in dataset box and update when changing config file
    def dataset_combobox_update(self):
        values = self.dataset_select['values']
        for value in values:
            if value == os.path.basename(os.path.normpath(self.raw_init_directory)):
               self.dataset_select.delete(0, tk.END)
               self.dataset_select.insert(0, value)
               
    # Method for changing datasets from the Configuration tab
    def update_dataset(self):
        dataset = self.dataset_select.get()
        self.remove_all_cuts()
        self.show_var.set(False)
        self.load_3d_bubble_data()
        if self.run_entry['state'] == tk.DISABLED:
            for child in self.bottom_frame_1.winfo_children():
                child.config(state=tk.NORMAL)

        try:
            split, end = dataset.split('-d', 1)
            values = self.ped_config_file_path_combobox['values']
            updated = False
            for value in values:
                if value.endswith(split + '-ped_config.txt'):
                    self.ped_config_file_path_combobox.delete(0, tk.END)
                    self.ped_config_file_path_combobox.insert(0, value)
                    updated = True
                    break
            if updated:
                self.reco_filename = 'reco_events.npy'
                self.config_combobox_update()
            ## if a config file is not found, this process will just set default values, error, and allow user to enter values
            else:
                self.piezo_combobox.delete(0, tk.END)
                self.piezo_combobox.insert(0, [self.piezo])
                self.piezo_combobox.current(0)
                self.dytran_combobox.delete(0, tk.END)
                self.dytran_combobox.insert(0, [self.dytran])
                self.dytran_combobox.current(0)
                self.piezo_selector_combobox.delete(0, tk.END)
                self.piezo_selector_combobox.insert(0, [self.piezo])
                self.piezo_selector_combobox.current(0)

                self.raw_init_directory = os.path.join(self.base_directory, dataset)
                self.log_directory = os.path.join(self.raw_init_directory, 'logs')
                self.raw_directory = self.raw_init_directory
                self.scan_directory = os.path.join(self.ped_directory, 'scan_output_' + dataset[:-5])
                self.npy_directory = os.path.join(self.ped_directory, 'npy', dataset[:-5])
                self.raw_init_directory_path_var.set(os.path.normpath(self.raw_init_directory))
                self.raw_directory = self.raw_init_directory
                self.scan_directory_path_var.set(self.scan_directory)
                self.npy_directory_path_var.set(self.npy_directory)

                self.do_handscan_checkbutton['state'] = tk.NORMAL

                self.ped_config_file_path_combobox.delete(0, tk.END)
                self.ped_config_file_path_var = os.path.join(self.config_file_directory, split + '-ped_config.txt')
                self.load_config_values(self.ped_config_file_path_var)

                self.raw_directory_config_entry.delete(0, tk.END)
                self.raw_directory_config_entry.insert(0, self.raw_init_directory)

                self.scratch_directory_config_entry.delete(0, tk.END)
                self.scratch_directory_config_entry.insert(0, os.path.join(self.ped_directory, 'scratch', getpass.getuser()))

                self.init_run_config_entry.delete(0, tk.END)
                self.init_run_config_entry.insert(0, self.init_run)

                self.plc_temp_config_entry.delete(0, tk.END)
                self.plc_temp_config_entry.insert(0, self.plc_temp_var)

                self.relative_path_to_images_config_entry.delete(0, tk.END)
                self.relative_path_to_images_config_entry.insert(0, self.images_relative_path)

                self.image_naming_convention_select.delete(0, tk.END)
                self.image_naming_convention_select.insert(0, self.image_naming_convention)

                self.num_cams_config_entry.delete(0, tk.END)
                self.num_cams_config_entry.insert(0, self.num_cams)

                self.image_orientation_select.delete(0, tk.END)
                self.image_orientation_select.insert(0, self.image_orientation)

                self.first_frame_config_entry.delete(0, tk.END)
                self.first_frame_config_entry.insert(0, self.first_frame)

                self.init_frame_config_entry.delete(0, tk.END)
                self.init_frame_config_entry.insert(0, self.init_frame)

                self.last_frame_config_entry.delete(0, tk.END)
                self.last_frame_config_entry.insert(0, self.last_frame)

                self.update_num_cams()
                self.load_reco()
                self.get_raw_events()
                self.reset_event()
                self.set_init_dataset()
                self.reset_images()
                self.load_3d_bubble_data()
                
        except FileNotFoundError:
            self.logger.error(
                'Raw events not found for this dataset. Please ensure that the raw_events.npy file is present')
            self.num_cams = 0
            self.update_num_cams()
            if self.load_fastDAQ_piezo_checkbutton_var.get():
                self.load_fastDAQ_piezo_checkbutton_var.set(False)
                self.load_fastDAQ_piezo()
            if self.load_initial_data_checkbutton_var.get():
                self.load_initial_data_checkbutton_var.set(False)
                self.load_fastDAQ_analysis()
            if self.load_dytran_checkbutton_var.get():
                self.load_dytran_checkbutton_var.set(False)
                self.load_fastDAQ_dytran()
            for child in self.bottom_frame_1.winfo_children():
                child.config(state=tk.DISABLED)
        
        self.reset_handscan()

    # Method for changing data directories
    def update_directories(self):
        try:
            self.raw_init_directory_path_var.set(os.path.normpath(self.raw_init_directory))
            self.raw_directory = self.raw_init_directory
            self.scan_directory_path_var.set(self.scan_directory)
            self.npy_directory_path_var.set(self.npy_directory)

            if self.npy_directory_label['state'] != tk.DISABLED:
                # if self.reco_version_combobox.get() == 'devel':
                #     self.npy_directory = self.npy_directory.replace('current', 'devel')
                # else:
                #     self.npy_directory = self.npy_directory.replace('devel', 'current')
                if self.reco_version_combobox.get() == 'abub':
                    self.npy_directory = os.path.join(self.ped_directory, 'npy', self.dataset, 'abub')
                elif self.reco_version_combobox.get() == 'devel':
                    self.npy_directory = os.path.join(self.ped_directory, 'npy', self.dataset, 'devel')
                else:
                    self.npy_directory = os.path.join(self.ped_directory, 'npy', self.dataset)
                self.npy_directory_path_var.set(self.npy_directory)
                self.load_reco()
                values = sorted(self.reco_events.dtype.names) if self.reco_events is not None else ('')
                self.add_display_var_combobox['values'] = values
                
            if (not os.path.exists(self.raw_directory)) or (not os.path.exists(self.npy_directory)):
                raise FileNotFoundError

            self.get_raw_events()
            self.load_reco()
            self.reset_event()
            self.remove_all_cuts()
            self.set_init_dataset()
            self.num_cams = int(self.num_cams_config_entry.get())
            self.update_num_cams()
            self.reset_images()

            if self.run_entry['state'] == tk.DISABLED:
                for child in self.bottom_frame_1.winfo_children():
                    child.config(state=tk.NORMAL)

        except FileNotFoundError:
            self.logger.error(
                'One or more directories not found at given paths.\nPlease check paths, and also ensure that reco_events.npy and raw_events.npy '
                'exist in their respective directories')
            self.num_cams = 0
            self.update_num_cams()
            if self.load_fastDAQ_piezo_checkbutton_var.get():
                self.load_fastDAQ_piezo_checkbutton_var.set(False)
                self.load_fastDAQ_piezo()
            if self.load_initial_data_checkbutton_var.get():
                self.load_initial_data_checkbutton_var.set(False)
                self.load_fastDAQ_analysis()
            if self.load_dytran_checkbutton_var.get():
                self.load_dytran_checkbutton_var.set(False)
                self.load_fastDAQ_dytran()
            for child in self.bottom_frame_1.winfo_children():
                child.config(state=tk.DISABLED)
        
        self.reset_handscan()

    # for when manual config path is updated
    def new_config_update(self):
        if os.path.exists(self.path_ped_config_directory.get()):
            self.config_file_directory = self.path_ped_config_directory.get()
            self.ped_config_file_path_combobox['values'] = self.get_configs()
        else:
            self.path_ped_config_directory.delete(0, tk.END)
            self.path_ped_config_directory.insert(0, self.config_file_directory)
            self.logger.error('Given config directory not found')

    # for when the config file path is changed
    def config_combobox_update(self):
        self.ped_config_file_path_var = os.path.join(self.config_file_directory, self.ped_config_file_path_combobox.get())
        self.remove_all_cuts()
        self.load_config_values(self.ped_config_file_path_var)

        self.piezo_combobox['values'] = [self.piezo]
        self.piezo_combobox.current(0)
        self.dytran_combobox['values'] = [self.dytran]
        self.dytran_combobox.current(0)
        self.piezo_selector_combobox['values'] = [self.piezo]
        self.piezo_selector_combobox.current(0)
        self.piezo_plot_t0_checkbutton_var.set(False)
        self.dytran_plot_t0_checkbutton_var.set(False)
        self.show_var.set(False)
        
        # update all of the widget values
        self.raw_init_directory_path_var.set(os.path.normpath(self.raw_init_directory))
        self.scan_directory_path_var.set(self.scan_directory)
        self.npy_directory_path_var.set(self.npy_directory)

        if self.npy_directory_label['state'] != tk.DISABLED:
            # if self.reco_version_combobox.get() == 'devel':
            #     self.npy_directory = self.npy_directory.replace('current', 'devel')
            # else:
            #     self.npy_directory = self.npy_directory.replace('devel', 'current')
            if self.reco_version_combobox.get() == 'abub':
                self.npy_directory = os.path.join(self.ped_directory, 'npy', self.dataset, 'abub')
                self.load_reco()
            elif self.reco_version_combobox.get() == 'devel':
                self.npy_directory = os.path.join(self.ped_directory, 'npy', self.dataset, 'devel')
                self.load_reco()
            else:
                self.npy_directory = os.path.join(self.ped_directory, 'npy', self.dataset)
                self.npy_directory_path_var.set(self.npy_directory)
                self.load_reco()
                values = sorted(self.reco_events.dtype.names) if self.reco_events is not None else ('')
                self.add_display_var_combobox['values'] = values
        
        self.raw_directory_config_entry.delete(0, tk.END)
        self.raw_directory_config_entry.insert(0, self.raw_init_directory)

        self.scratch_directory_config_entry.delete(0, tk.END)
        self.scratch_directory_config_entry.insert(0, self.extraction_path)

        self.init_run_config_entry.delete(0, tk.END)
        self.init_run_config_entry.insert(0, self.init_run)

        self.plc_temp_config_entry.delete(0, tk.END)
        self.plc_temp_config_entry.insert(0, self.plc_temp_var)

        self.relative_path_to_images_config_entry.delete(0, tk.END)
        self.relative_path_to_images_config_entry.insert(0, self.images_relative_path)

        self.image_naming_convention_select.delete(0, tk.END)
        self.image_naming_convention_select.insert(0, self.image_naming_convention)

        self.num_cams_config_entry.delete(0, tk.END)
        self.num_cams_config_entry.insert(0, self.num_cams)

        self.image_orientation_select.delete(0, tk.END)
        self.image_orientation_select.insert(0, self.image_orientation)

        self.first_frame_config_entry.delete(0, tk.END)
        self.first_frame_config_entry.insert(0, self.first_frame)

        self.init_frame_config_entry.delete(0, tk.END)
        self.init_frame_config_entry.insert(0, self.init_frame)

        self.last_frame_config_entry.delete(0, tk.END)
        self.last_frame_config_entry.insert(0, self.last_frame)

        self.piezo_config_entry.delete(0, tk.END)
        self.piezo_config_entry.insert(0, self.piezo)

        self.dytran_config_entry.delete(0, tk.END)
        self.dytran_config_entry.insert(0, self.dytran)

        self.jar_radius_config_entry.delete(0, tk.END)
        self.jar_radius_config_entry.insert(0, self.radius)

        self.jar_positive_z_config_entry.delete(0, tk.END)
        self.jar_positive_z_config_entry.insert(0, self.positive_z)

        self.jar_negative_z_config_entry.delete(0, tk.END)
        self.jar_negative_z_config_entry.insert(0, self.negative_z)

        self.update_num_cams()
        self.get_raw_events()

        self.reset_event()
        self.remove_all_cuts()
        self.reset_cuts()
        self.set_init_dataset()
        self.load_3d_bubble_data()
        self.reset_images()
        self.dataset_combobox_update()
        
        self.reset_handscan()

        values = sorted(self.reco_events.dtype.names) if self.reco_events is not None else ('')
        self.add_display_var_combobox['values'] = values

    # When saving a new config file, checks that all changes have been applied
    def are_unapplied_config_changes(self):
        if self.raw_init_directory != self.raw_directory_config_entry.get():
            return True
        if os.path.basename(self.scratch_directory_config_entry.get()) != getpass.getuser():
            return True
        if self.init_run != self.init_run_config_entry.get():
            return True
        if self.plc_temp_var != self.plc_temp_config_entry.get():
            return True
        if self.images_relative_path != self.relative_path_to_images_config_entry.get():
            return True
        if self.image_naming_convention != self.image_naming_convention_select.get():
            return True
        if self.num_cams != int(self.num_cams_config_entry.get()):
            return True
        if self.image_orientation != self.image_orientation_select.get():
            return True
        if self.init_frame != self.init_frame_config_entry.get():
            return True
        if self.first_frame != self.first_frame_config_entry.get():
            return True
        if self.last_frame != self.last_frame_config_entry.get():
            return True
        if self.piezo != self.piezo_config_entry.get():
            return True
        if self.dytran != self.dytran_config_entry.get():
            return True
        if self.radius != self.jar_radius_config_entry.get():
            return True
        if self.positive_z != self.jar_positive_z_config_entry.get():
            return True
        if self.negative_z != self.jar_negative_z_config_entry.get():
            return True
        
        return False
    
    # for when values are changed without updating config paths
    def update_vars_config(self):
        #raw_init_directory, sractch_directory, init_run
        if os.path.basename(self.scratch_directory_config_entry.get()) == getpass.getuser():
            self.extraction_path = self.scratch_directory_config_entry.get()
        else:
            self.scratch_directory_config_entry.delete(0, tk.END)
            self.scratch_directory_config_entry.insert(0, self.extraction_path)
            self.logger.error("Error: Scratch Directory must be a Folder matching the current User Name. \n\nEx:C:\\Users\\User_Name\\PICO\\EventDisplay\\scratch\\User_Name\n\nPlease change Scratch Directory to a Folder matching the current User Name\n\nNo changes applied.  ")
        self.raw_init_directory = self.raw_directory_config_entry.get()
        self.init_run = self.init_run_config_entry.get()

        # plc var
        self.plc_temp_var = self.plc_temp_config_entry.get()
        self.load_plc_text()

        # image related vars
        self.image_orientation = self.image_orientation_select.get()
        self.num_cams = int(self.num_cams_config_entry.get())
        self.images_relative_path = self.relative_path_to_images_config_entry.get()
        self.image_naming_convention = self.image_naming_convention_select.get()
        self.init_frame = self.init_frame_config_entry.get()
        self.first_frame = self.first_frame_config_entry.get()
        self.last_frame = self.last_frame_config_entry.get()
        self.piezo = self.piezo_config_entry.get()
        self.dytran = self.dytran_config_entry.get()
        self.radius = self.jar_radius_config_entry.get()
        self.positive_z = self.jar_positive_z_config_entry.get()
        self.negative_z = self.jar_negative_z_config_entry.get()

        # reset piezo & dytran combo boxes
        self.piezo_combobox['values'] = [self.piezo]
        self.piezo_combobox.current(0)
        self.dytran_combobox['values'] = [self.dytran]
        self.dytran_combobox.current(0)
        self.piezo_selector_combobox['values'] = [self.piezo]
        self.piezo_selector_combobox.current(0)
        
        # reset other values
        self.show_all_reco_var.set(False)
        self.load_3d_bubble_data()
        self.update_num_cams()
        self.reset_images()
        self.remove_all_cuts()
        self.reset_event()

    def update_num_cams(self):
        # reset the number of canvases
        for canvas in self.canvases:
            canvas.delete('all')
        self.canvases = []
        for cam in range(0, self.num_cams):
            canvas = tk.Canvas(self.camera_tab, width=self.init_image_width, height=self.init_image_height)
            canvas.bind('<ButtonPress-1>', self.on_button_press)
            canvas.zoom = 0

            canvas.image = canvas.create_image(0, 0, anchor=tk.NW, image=None)
            canvas.bottom_text = canvas.create_text(10, self.init_image_height - 20, anchor=tk.NW, text='',
                                                    fill='red')
            canvas.grid(row=0, column=1 * cam, columnspan=1, sticky='NW')
            canvas.cam = cam
            self.canvases.append(canvas)

    def toggle_reco_widgets(self, state):
        self.draw_crosshairs_button.config(state=state)
        self.plot_t0_checkbutton.config(state=state)
        self.piezo_plot_t0_checkbutton.config(state=state)
        self.dytran_plot_t0_checkbutton.config(state=state)
        for child in self.bottom_frame_2.winfo_children():
            child.config(state=state)
        self.convert_reco_button['state'] = tk.NORMAL
        self.select_reco_button['state'] = tk.NORMAL
        self.reco_label['state'] = tk.NORMAL
        self.reco_availability_label['state'] = tk.NORMAL

    # Creates a new config file based on combo-box/text inputs that saves with a record of the user, day, and time
    def save_current_config(self):
        if self.are_unapplied_config_changes():
            messagebox.showwarning('Alert','You have unapplied changes to the config vars. Applying before saving.')
            
        self.update_vars_config()       
        if self.dataset_select.get().split('-d', 1)[0] not in self.ped_config_file_path_combobox.get(): # This compares values in comboboxes for Dataset and path to config file. 
           messagebox.showwarning('Alert','You have changed dataset or config file selection without updating. No file has been created.')
           return
        elif self.dataset_select.get().split('-d', 1)[0] not in self.ped_config_file_path_var: #This compares value in the dataset Combobox with the currently selected config file, that way we do not create a config file with title 30l-16 with data from 2l-15 accidently for example.
           messagebox.showwarning('Alert','You have changed dataset or config file selection without updating. No file has been created.')
           return
        new_config_file_name = '{}-ped_config_{}_{}.txt'.format(self.dataset_select.get().split('-d', 1)[0], getpass.getuser(), time.strftime('%a_%b_%d_%H_%M_%S_%Y'))
        self.new_config_file_path = os.path.join(self.config_file_directory, new_config_file_name)
        self.new_file_creator = open(self.new_config_file_path, 'w+')
        if self.image_naming_convention == 'cam0_image0.png':
           image_convention = 0
        else:
           image_convention = 1
        if self.image_orientation == '0':
           orientation = 0
        elif self.image_orientation == '90':
           orientation = 1
        elif self.image_orientation == '180':
           orientation = 2
        else:
           orientation = 3  
        text_in_file = "Raw Directory:\n{}\nScratch Directory:\n{}\nInitial Run:\n{}\nPLC temperature var:\n{}\nRelative Path to images:\n{}\nImage naming convention(0 for 'cam0_image0.png', 1 for 'cam0image  0.bmp'):\n{}\nNumber of Cameras:\n{}\nImage Orientation (0 for '0', 1 for '90', 2 for '180', 3 for '270'):\n{}\nfirst frame:\n{}\ntrig frame:\n{}\nlast frame:\n{}\nPiezo:\n{}\nDytran:\n{}\nJar - Radius (mm):\n{}\nJar - Positive Z Height (mm):\n{}\nJar - Negative Z Height (mm):\n{}".format(self.raw_init_directory, os.path.dirname(self.extraction_path), self.init_run, self.plc_temp_var, self.images_relative_path, image_convention, self.num_cams, orientation, self.first_frame, self.init_frame, self.last_frame, self.piezo, self.dytran, self.radius, self.positive_z, self.negative_z)
        self.new_file_creator.write(text_in_file)
        self.new_file_creator.close()
        self.new_config_update()

    # resets handscan-widget after loading/applying new config file details
    def reset_handscan(self):
        if self.do_handscan_checkbutton_var.get():
            messagebox.showwarning('Alert','You have changed dataset while handscanning')

        self.do_handscan_checkbutton['state'] = tk.NORMAL
        self.do_handscan_checkbutton_var.set(0)
        flag=tk.DISABLED
        for i in range(0, 8):
            self.nbub_button[i].config(state=flag)

        for button in self.source_buttons:
            button.config(state=flag)

        self.isgoodtrigger_button.config(state=flag)
        self.crosshairsgood_button.config(state=flag)
        self.comment_label.config(state=flag)
        self.comment_entry.config(state=flag)
        self.submit_scan_button.config(state=flag)

    def use_custom_reco(self):
        if self.use_custom_reco_checkbutton_var.get():
            self.select_custom_reco_combobox['state'] = tk.NORMAL
            self.select_custom_reco_combobox['values'] = self.get_custom_reco()
            self.update_custom_reco_button['state'] = tk.NORMAL
        else:
            self.select_custom_reco_combobox['state'] = tk.DISABLED
            self.update_custom_reco_button['state'] = tk.DISABLED

    def update_custom_reco(self):
        reco_filename = self.select_custom_reco_combobox.get()

        if os.path.isfile(reco_filename):
            reco_directory = os.path.split(reco_filename)[0]
            if not reco_directory == self.npy_directory:
                directory_response = messagebox.askyesno('Selected File Not in NPY Directory', 'Current Dataset is: {}. \nSelected reco file is not in the NPY Directory for the current dataset. \nWould you like still like to use the selected reco_events file?'.format(self.dataset))
                if directory_response == 1:
                    self.reco_filename = os.path.split(reco_filename)[1]
                else:
                    return
            else:
                self.reco_filename = os.path.split(reco_filename)[1]
        else:
            return

        self.load_reco()
        self.load_3d_bubble_data()
        self.update_num_cams()
        self.reset_images()
        self.remove_all_cuts()
        self.reset_event()

    #Set-up frames, boxes, and labels for the tab
    def create_configuration_widgets(self):
        self.config_tab = tk.Frame(self.notebook)
        self.notebook.add(self.config_tab, text='Configuration')

        # configuration tab
        # setup frames within tab
        self.config_tab_main = tk.Frame(self.config_tab, bd=5, relief=tk.SUNKEN)
        self.config_tab_main.grid(row=0, column=0, sticky='NW')

        self.config_tab_vars = tk.Frame(self.config_tab, bd=5, relief=tk.SUNKEN)
        self.config_tab_vars.grid(row=4, column=0, sticky='NW')

        # frame 1
        self.dataset_label = tk.Label(self.config_tab_main, text='Dataset:', width=20)
        self.dataset_label.grid(row=0, column=0, sticky='WE')
        self.dataset_select = ttk.Combobox(self.config_tab_main, values=self.get_datasets(), width=122)
        self.dataset_select.grid(row=0, column=1, sticky='WE')
        self.set_init_dataset()

        self.update_dataset_button = tk.Button(self.config_tab_main, text='Update Dataset', command=self.update_dataset)
        self.update_dataset_button.grid(row=0, column=2, sticky='NW')

        self.raw_directory_label = tk.Label(self.config_tab_main, text='Raw Directory:')
        self.raw_directory_label.grid(row=1, column=0, sticky='WE')
        self.raw_directory_path = tk.Label(self.config_tab_main, textvariable=self.raw_init_directory_path_var)
        self.raw_directory_path.grid(row=1, column=1, sticky='WE')

        self.scan_directory_label = tk.Label(self.config_tab_main, text='Scan Directory:')
        self.scan_directory_label.grid(row=2, column=0, sticky='WE')
        self.scan_directory_path = tk.Label(self.config_tab_main, textvariable=self.scan_directory_path_var)
        self.scan_directory_path.grid(row=2, column=1, sticky='WE')

        self.npy_directory_label = tk.Label(self.config_tab_main, text='.npy Directory:')
        self.npy_directory_label.grid(row=3, column=0, sticky='WE')
        self.npy_directory_path = tk.Label(self.config_tab_main, textvariable=self.npy_directory_path_var)
        self.npy_directory_path.grid(row=3, column=1, sticky='WE')

        self.reco_version_label = tk.Label(self.config_tab_main, text='Reco Version:')
        self.reco_version_label.grid(row=4, column=0, sticky='WE')
        self.reco_version_combobox = ttk.Combobox(self.config_tab_main, values=['current', 'devel', 'abub'])
        self.reco_version_combobox.grid(row=4, column=1, sticky='WE')

        self.update_directory_button = tk.Button(self.config_tab_main, text='Apply',
                                                 command=self.update_directories)
        self.update_directory_button.grid(row=4, column=2, sticky='WE')

        self.use_custom_reco_checkbutton = tk.Checkbutton(self.config_tab_main,
                                                       text='Use Custom Reco',
                                                       variable=self.use_custom_reco_checkbutton_var,
                                                       command=self.use_custom_reco)
        self.use_custom_reco_checkbutton.grid(row=5, column=0, sticky='WE')
        self.select_custom_reco_combobox  = ttk.Combobox(self.config_tab_main, values=self.get_custom_reco)
        self.select_custom_reco_combobox['state'] = tk.DISABLED
        self.select_custom_reco_combobox.grid(row=5, column=1, sticky='WE')

        self.update_custom_reco_button = tk.Button(self.config_tab_main, text='Apply Custom Reco',
                                                 command=self.update_custom_reco)
        self.update_custom_reco_button['state'] = tk.DISABLED
        self.update_custom_reco_button.grid(row=5, column=2, sticky='WE')

        # frame 2
        self.path_ped_config_directory_label = tk.Label(self.config_tab_vars, text='Path to config directory')
        self.path_ped_config_directory_label.grid(row=0, column=0, sticky='WE')
        self.path_ped_config_directory = tk.Entry(self.config_tab_vars, width=125)
        self.path_ped_config_directory.grid(row=0, column=1, columnspan=5, sticky='NW')

        self.update_vars_config_button = tk.Button(self.config_tab_vars, text='Update Config Directory',
                                                   command=self.new_config_update)
        self.update_vars_config_button.grid(row=0, column=6, sticky='NW')

        self.update_vars_config_button = tk.Button(self.config_tab_vars, text='Apply Current Values',
                                                   command=self.update_vars_config)
        self.update_vars_config_button.grid(row=2, column=6, sticky='NW')

        self.ped_config_file_path = tk.Label(self.config_tab_vars, text='Config.txt file:')
        self.ped_config_file_path.grid(row=1, column=0, sticky='WE')
        self.ped_config_file_path_combobox = ttk.Combobox(self.config_tab_vars, width=60)
        self.ped_config_file_path_combobox.grid(row=1, column=1, columnspan=5, sticky='WE')

        self.update_config_combobox_button = tk.Button(self.config_tab_vars, text='Select Config File',
                                                       command=self.config_combobox_update)
        self.update_config_combobox_button.grid(row=1, column=6, sticky='NW')
        
        self.save_config_button = tk.Button(self.config_tab_vars, text='Save As Config File', bg = "#CD5C5C",
                                                       command=self.save_current_config)
        self.save_config_button.grid(row=3, column=6, sticky='NW')

        self.raw_directory_config_label = tk.Label(self.config_tab_vars, text='Raw Directory:')
        self.raw_directory_config_label.grid(row=2, column=0, sticky='WE')
        self.raw_directory_config_entry = tk.Entry(self.config_tab_vars, width=30)
        self.raw_directory_config_entry.insert(0, self.raw_init_directory)
        self.raw_directory_config_entry.grid(row=2, column=1, columnspan=5, sticky='WE')

        self.scratch_directory_config_label = tk.Label(self.config_tab_vars, text='Scratch Directory:')
        self.scratch_directory_config_label.grid(row=3, column=0, sticky='WE')
        self.scratch_directory_config_entry = tk.Entry(self.config_tab_vars, width=30)
        self.scratch_directory_config_entry.insert(0, self.extraction_path)
        self.scratch_directory_config_entry.grid(row=3, column=1, columnspan=5, sticky='WE')

        self.init_run_config_label = tk.Label(self.config_tab_vars, text='Initial Run:')
        self.init_run_config_label.grid(row=4, column=0, sticky='WE')
        self.init_run_config_entry = tk.Entry(self.config_tab_vars, width=10)
        self.init_run_config_entry.insert(0, self.init_run)
        self.init_run_config_entry.grid(row=4, column=1, sticky='WE')

        self.plc_temp_config_label = tk.Label(self.config_tab_vars, text='PLC temperature var:')
        self.plc_temp_config_label.grid(row=4, column=2, sticky='WE')
        self.plc_temp_config_entry = tk.Entry(self.config_tab_vars, width=10)
        self.plc_temp_config_entry.insert(0, self.plc_temp_var)
        self.plc_temp_config_entry.grid(row=4, column=3,  sticky='WE')

        self.relative_path_to_images_config_label = tk.Label(self.config_tab_vars, text='Relative path to images:')
        self.relative_path_to_images_config_label.grid(row=4, column=4, sticky='WE')
        self.relative_path_to_images_config_entry = tk.Entry(self.config_tab_vars, width=10)
        self.relative_path_to_images_config_entry.insert(0, self.images_relative_path)
        self.relative_path_to_images_config_entry.grid(row=4, column=5, sticky='WE')

        self.image_naming_convention_label = tk.Label(self.config_tab_vars, text='Image naming convention:')
        self.image_naming_convention_label.grid(row=5, column=0, sticky='WE')
        self.image_naming_convention_select = ttk.Combobox(self.config_tab_vars, values=self.image_naming_conventions)
        self.image_naming_convention_select.insert(0, self.image_naming_convention)
        self.image_naming_convention_select.grid(row=5, column=1, sticky='WE')

        self.num_cams_config_label = tk.Label(self.config_tab_vars, text='Number of cameras:')
        self.num_cams_config_label.grid(row=5, column=2, sticky='WE')
        self.num_cams_config_entry = tk.Entry(self.config_tab_vars, width=10)
        self.num_cams_config_entry.insert(0, self.num_cams)
        self.num_cams_config_entry.grid(row=5, column=3, sticky='WE')

        self.image_orientation_config_label = tk.Label(self.config_tab_vars, text='Image orientation:')
        self.image_orientation_config_label.grid(row=5, column=4, sticky='WE')
        self.image_orientation_select = ttk.Combobox(self.config_tab_vars, values=self.image_orientations)
        self.image_orientation_select.insert(0, self.image_orientation)
        self.image_orientation_select.grid(row=5, column=5, sticky='WE')

        self.first_frame_config_label = tk.Label(self.config_tab_vars, text='first frame:')
        self.first_frame_config_label.grid(row=6, column=0, sticky='WE')
        self.first_frame_config_entry = tk.Entry(self.config_tab_vars, width=10)
        self.first_frame_config_entry.insert(0, self.first_frame)
        self.first_frame_config_entry.grid(row=6, column=1, sticky='WE')

        self.init_frame_config_label = tk.Label(self.config_tab_vars, text='trig frame:')
        self.init_frame_config_label.grid(row=6, column=2, sticky='WE')
        self.init_frame_config_entry = tk.Entry(self.config_tab_vars, width=10)
        self.init_frame_config_entry.insert(0, self.init_frame)
        self.init_frame_config_entry.grid(row=6, column=3, sticky='WE')

        self.last_frame_config_label = tk.Label(self.config_tab_vars, text='last frame:')
        self.last_frame_config_label.grid(row=6, column=4, sticky='WE')
        self.last_frame_config_entry = tk.Entry(self.config_tab_vars, width=10)
        self.last_frame_config_entry.insert(0, self.last_frame)
        self.last_frame_config_entry.grid(row=6, column=5, sticky='WE')

        self.piezo_config_label = tk.Label(self.config_tab_vars, text='Piezo:')
        self.piezo_config_label.grid(row=7, column=0, sticky='WE')
        self.piezo_config_entry = tk.Entry(self.config_tab_vars, width=10)
        self.piezo_config_entry.insert(0, self.piezo)
        self.piezo_config_entry.grid(row=7, column=1, sticky='WE')

        self.dytran_config_label = tk.Label(self.config_tab_vars, text='Dytran:')
        self.dytran_config_label.grid(row=7, column=2, sticky='WE')
        self.dytran_config_entry = tk.Entry(self.config_tab_vars, width=10)
        self.dytran_config_entry.insert(0, self.dytran)
        self.dytran_config_entry.grid(row=7, column=3, sticky='WE')

        self.jar_radius_config_label = tk.Label(self.config_tab_vars, text='Jar - Radius (mm):')
        self.jar_radius_config_label.grid(row=8, column=0, sticky='WE')
        self.jar_radius_config_entry = tk.Entry(self.config_tab_vars, width=10)
        self.jar_radius_config_entry.insert(0, self.radius)
        self.jar_radius_config_entry.grid(row=8, column=1, sticky='WE')

        self.jar_positive_z_config_label = tk.Label(self.config_tab_vars, text='Jar - Positive Z Height (mm):')
        self.jar_positive_z_config_label.grid(row=8, column=2, sticky='WE')
        self.jar_positive_z_config_entry = tk.Entry(self.config_tab_vars, width=10)
        self.jar_positive_z_config_entry.insert(0, self.positive_z)
        self.jar_positive_z_config_entry.grid(row=8, column=3, sticky='WE')

        self.jar_negative_z_config_label = tk.Label(self.config_tab_vars, text='Jar - Negative Z Height (mm):')
        self.jar_negative_z_config_label.grid(row=8, column=4, sticky='WE')
        self.jar_negative_z_config_entry = tk.Entry(self.config_tab_vars, width=10)
        self.jar_negative_z_config_entry.insert(0, self.negative_z)
        self.jar_negative_z_config_entry.grid(row=8, column=5, sticky='WE')


