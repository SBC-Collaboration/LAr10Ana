# Imports
import subprocess, os, platform
import tkinter as tk
from PIL import Image, ImageChops, ImageOps, ImageTk
import cv2


class Camera(tk.Frame):
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)

        # Initial Functions
        self.create_camera_widgets()
        self.zip_flag

    def reset_images(self):
        self.load_event_sbc()
        self.load_plc_text()
        self.load_fastDAQ_piezo()
        self.load_slowDAQ()
        # self.load_fastDAQ_dytran()
        self.load_fastdaq_scintillation()
        self.load_fastDAQ_analysis()
        self.frame = self.init_frame
        self.diff_checkbutton_var.set(False)
        self.invert_checkbutton_var.set(False)

        for canvas in self.canvases:
            self.reset_zoom(canvas)

        self.update_images()

    def get_image_path(self, cam, frame):
        print(self.image_directory)
        if self.image_naming_convention == self.image_naming_conventions[0]:
            path = os.path.join(self.image_directory, 'cam{}_image{}.png'.format(cam, frame))
        elif self.image_naming_convention == self.image_naming_conventions[1]:
            # handle the leading spaces in the image names
            frame = '{:>3}'.format(frame)
            path = os.path.join(self.image_directory, 'cam{}image{}.bmp'.format(cam, frame))
        elif self.image_naming_convention == self.image_naming_conventions[2]:
            # handle the leading zeros in the image names
            frame = str(frame).zfill(2)
            # camera numbering starts at 1, so cam + 1
            path = os.path.join(self.image_directory, 'cam{}-img{}.png'.format(cam + 1, frame))
        else:
            path = os.path.join(self.image_directory, 'cam{}_image{}.png'.format(cam, frame))
            self.error += ('Image naming convention not found\n')

        return path

    def update_images(self):
        error = ' '
        for canvas in self.canvases:
            path = self.get_image_path(canvas.cam, self.frame)
            image = self.load_image(path, canvas)

            zoom = '{:.1f}'.format(canvas.image_width / self.native_image_width)
            if self.diff_checkbutton_var.get():
                path = self.get_image_path(canvas.cam, self.first_frame)
                first_frame = self.load_image(path, canvas)

                image = ImageOps.autocontrast(ImageChops.difference(first_frame, image))
                
                template = 'frame: {} zoom: {}x (diff wrt {})                  {}/{}'
                bottom_text = template.format(self.frame, zoom, self.first_frame, self.run, self.event)
            else:
                template = 'frame: {} zoom: {}x                                   {}/{}'
                bottom_text = template.format(self.frame, zoom, self.run, self.event)

            if self.invert_checkbutton_var.get():
                image = ImageOps.invert(image)

            canvas.photo = ImageTk.PhotoImage(image)
            canvas.itemconfig(canvas.image, image=canvas.photo)
            canvas.itemconfig(canvas.bottom_text, text=bottom_text)

        self.draw_crosshairs()

    def load_image(self, path, canvas):
        try:
            path = self.zipped_event.open(path) if self.zip_flag else path
            image = Image.open(path)
            if self.image_orientation == '90':
                image = image.transpose(Image.ROTATE_90)
            if self.image_orientation == '180':
                image = image.transpose(Image.ROTATE_180)
            if self.image_orientation == '270':
                image = image.transpose(Image.ROTATE_270)
        except (KeyError, FileNotFoundError):
            self.logger.info('Did not find image at {}'.format(path))
            image = Image.open(os.path.join(self.ped_directory, 'notfound.jpeg'))
        except IOError:
            self.error += ('image format problem, attempting to recover\n')
            cv2img = cv2.imread(path)
            cv2.imwrite('ped_temp.jpg', cv2img)
            image = Image.open('ped_temp.jpg')

        # Image in RGBA, Pillow cannot deal with that, so convert to RGB
        # print('Image mode: ', image.mode)
        if image.mode == 'RGBA':
            r, g, b, a = image.split()
            image = Image.merge('RGB', (r, g, b))
        if image.mode == 'P':
            image = image.convert('L')

        self.native_image_width, self.native_image_height = image.size
        image = image.resize((int(canvas.image_width), int(canvas.image_height)),
                             self.antialias_checkbutton_var.get())
        image = image.crop((canvas.crop_left, canvas.crop_bottom, canvas.crop_right, canvas.crop_top))

        return image

    def reset_zoom(self, canvas):
        canvas.zoom = 0
        canvas.crop_left = 0
        canvas.crop_bottom = 0
        canvas.crop_right = self.init_image_width
        canvas.crop_top = self.init_image_height
        canvas.image_width = self.init_image_width
        canvas.image_height = self.init_image_height

    def on_button_press(self, event):
        canvas = event.widget
        old_width = canvas.image_width
        old_height = canvas.image_height

        canvas.zoom += 1
        if canvas.zoom > self.max_zoom:
            self.reset_zoom(canvas)
        else:
            canvas.image_width = 2 ** (canvas.zoom - 1) * self.native_image_width
            canvas.image_height = 2 ** (canvas.zoom - 1) * self.native_image_height
            if self.native_image_width < self.init_image_width:
                canvas.image_width = 2 ** (canvas.zoom) * self.init_image_width
            if self.native_image_height < self.init_image_height:
                canvas.image_height = 2 ** (canvas.zoom) * self.init_image_height

            new_center_x = (event.x + canvas.crop_left) * (canvas.image_width / old_width)
            new_center_y = (event.y + canvas.crop_bottom) * (canvas.image_height / old_height)
            if new_center_x < self.init_image_width / 2:
                # click was too far left, not enough new image for center to be here
                new_center_x = self.init_image_width / 2
            if new_center_x + self.init_image_width / 2 > canvas.image_width:
                # click was too far right, not enough new image for center to be here
                new_center_x = canvas.image_width - self.init_image_width / 2
            if new_center_y < self.init_image_height / 2:
                # click was too far up, not enough new image for center to be here
                new_center_y = self.init_image_height / 2
            if new_center_y + self.init_image_height / 2 > canvas.image_height:
                # click was too far down, not enough new image for center to be here
                new_center_y = canvas.image_height - self.init_image_height / 2

            canvas.crop_left = new_center_x - self.init_image_width / 2
            canvas.crop_bottom = new_center_y - self.init_image_height / 2
            canvas.crop_right = new_center_x + self.init_image_width / 2
            canvas.crop_top = new_center_y + self.init_image_height / 2

        self.update_images()

    def make_video(self):
        if self.zip_flag:
            self.logger.error('Currently navigating unextracted run in zip format. This feature is only available if the run has been extracted either in ordinary raw data directory or in users scratch space')            
            return

        if self.image_naming_convention == self.image_naming_conventions[0]:
            img_conv = '_image%d.png'
        elif self.image_naming_convention == self.image_naming_conventions[1]:
            img_conv = 'image  %d.bmp'
        else:
            print("ERROR: unknown image naming convention")
            return
        
        camnum = int(self.make_video_entry.get())
        tmp_folder_path = os.path.join( self.extraction_path, 'tmp')
        tmp_file_path = os.path.join( tmp_folder_path, "tmp_cam" + str(camnum) + '_' + self.dataset + '_' + self.run + '_' + str(self.event) + ".mp4")
        out_file_path = os.path.join( tmp_folder_path, "cam"     + str(camnum) + '_' + self.dataset + '_' + self.run + '_' + str(self.event) + ".mp4")
        run_folder_path = os.path.join(self.raw_directory, self.run)
        event_folder_path = os.path.join(run_folder_path, str(self.event))
        event_folder_path = os.path.join(event_folder_path, self.images_relative_path )
        img_path = os.path.join( event_folder_path, 'cam' + str(camnum) + img_conv )
        img_path = '\"' + img_path + '\"'
        # print( "make video image path: ", img_path )
        # print( "image orientation: ", self.image_orientation)
        # print( str(self.first_frame) )

        # if os.path.exists(tmp_file_path):
        #     os.system('rm ' + tmp_file_path)
        if not os.path.exists(out_file_path):
            # os.system('rm ' + out_file_path)
            try:
                os.system( "ffmpeg -start_number " + str(self.first_frame) + " -r 1/0.2 -i " + img_path + " -c:v libx264 -r 30 -pix_fmt yuv420p " + tmp_file_path)
                if self.image_orientation == '270':
                    os.system( "ffmpeg -i " + tmp_file_path + " -vf \"transpose=1\" -s 400x624 " + out_file_path )
                elif self.image_orientation == '180':
                    os.system( "ffmpeg -i " + tmp_file_path + " -vf \"transpose=1,transpose=1\" -s 400x624 " + out_file_path )
                elif self.image_orientation == '90':
                    os.system( "ffmpeg -i " + tmp_file_path + " -vf \"transpose=2\" -s 400x624 " + out_file_path )
                else:
                    os.system( "ffmpeg -i " + tmp_file_path + " -s 400x624 -c:a copy " + out_file_path )
            except:
                print("ffmpeg ERROR -- failed to make video")

        try:
            # os.system("open " + out_file_path)
            if platform.system() == 'Darwin':       # macOS
                subprocess.call(('open', out_file_path))
            elif platform.system() == 'Windows':    # Windows
                os.startfile(out_file_path)
            else:                                   # linux variants
                # subprocess.call(('xdg-open', out_file_path))
                # subprocess.call(('cvlc', '--no-audio', out_file_path))
                subprocess.call(('ffplay', '-loop', '0', out_file_path))
        except:
            print("Video player ERROR -- failed to open video")
    
    def draw_crosshairs(self):
        for canvas in self.canvases:
            canvas.delete('crosshair')

        if not self.draw_crosshairs_var.get() or not self.reco_row:  # no reco row means we don't have reco data
            return

        try:
            if self.reco_row['nbub'] < 1:
                return
        except:
            print('WARNING: no nbub var')
            return

        for ibub in range(1, self.reco_row['nbub'] + 1):
            # print('crosshair self.run pre-recorow: ', self.run)
            # print('crosshair reco.run pre-recorow: ', self.reco_row['run'])
            # print('crosshair self.ev pre-recorow: ', self.event)
            # print('crosshair reco.ev pre-recorow: ', self.reco_row['ev'])
            # print('crosshair reco.nbub pre-recorow: ', self.reco_row['nbub'])
            self.load_reco_row(ibub)
            # print('  crosshair self.run post-recorow: ', self.run)
            # print('  crosshair reco.run post-recorow: ', self.reco_row['run'])
            # print('  crosshair self.ev post-recorow: ', self.event)
            # print('  crosshair reco.ev post-recorow: ', self.reco_row['ev'])
            # print('  crosshair reco.nbub post-recorow: ', self.reco_row['nbub'])
            for canvas in self.canvases:
                x_zoom = canvas.image_width / self.native_image_width
                y_zoom = canvas.image_height / self.native_image_height

                bubble_x = self.reco_row['hori{}'.format(canvas.cam)]
                bubble_y = self.reco_row['vert{}'.format(canvas.cam)]

                x = canvas.image_width - (bubble_x + canvas.crop_left / x_zoom) * x_zoom
                y = (bubble_y - canvas.crop_bottom / y_zoom) * y_zoom

                if self.image_orientation == '0':
                    x = (bubble_x - canvas.crop_left / x_zoom) * x_zoom
                    y = canvas.image_height - (bubble_y + canvas.crop_bottom / y_zoom) * y_zoom
                
                # print(' nbub: ', self.reco_row['nbub'])
                # print(' ibub: ', ibub)
                # print(' crosshair coord: ', str(x), str(y))
                
                canvas.create_line(x - 11, y, x - 5, y, fill='red', tag='crosshair')
                canvas.create_line(x + 5, y, x + 11, y, fill='red', tag='crosshair')
                canvas.create_line(x, y - 11, x, y - 5, fill='red', tag='crosshair')
                canvas.create_line(x, y + 5, x, y + 11, fill='red', tag='crosshair')
                canvas.create_oval(x - 8, y - 8, x + 8, y + 8, outline='red', tag='crosshair')

    def change_nbub(self):
        if self.nbub_button_var.get() > 1:
            for button in self.source_buttons:
                button.config(state=tk.ACTIVE)
            # self.source_button_var.set(0)

    # Moving forwards and backwards through image frames
    def load_frame(self, frame):
        self.frame = str(frame)

        path = self.get_image_path(0, self.frame)
        if self.zip_flag:
            try:
                self.zipped_event.open(path)
            except:
                self.frame = self.init_frame
        elif not os.path.isfile(path):
            self.frame = self.init_frame

        self.update_images()

    def create_camera_widgets(self):
        self.camera_tab = tk.Frame(self.notebook)
        self.notebook.add(self.camera_tab, text='Camera')

        # Cameras tab
        self.canvases = []
        for cam in range(0, self.num_cams):
            canvas = tk.Canvas(self.camera_tab, width=self.init_image_width, height=self.init_image_height)
            canvas.bind('<ButtonPress-1>', self.on_button_press)
            canvas.zoom = 0

            canvas.image = canvas.create_image(0, 0, anchor=tk.NW, image=None)
            canvas.bottom_text = canvas.create_text(10, self.init_image_height - 20, anchor=tk.NW, text='', fill='red')
            canvas.grid(row=0, column=1 * cam, columnspan=1, sticky='NW')
            canvas.cam = cam
            self.canvases.append(canvas)
