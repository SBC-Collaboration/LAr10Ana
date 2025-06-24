# Imports
import os
import re
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import pandas as pd
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from glob import glob
#
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from scipy.optimize import curve_fit
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import math


class ThreeDBubble(tk.Frame):
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)

        self.show_all_reco_var = tk.BooleanVar(value=False)
        self.show_var = tk.BooleanVar(value=False)

        self.X_var = 0
        self.Y_var = 0
        self.Z_var = 0
        self.edge_var = 0

        # Initial Functions
        self.create_three_d_bubble_widgets()
        self.three_d_bubble_canvas_setup()

    def three_d_bubble_canvas_setup(self):

        # Create Figures, Axes, and Canvases
        self.three_d_bubble_fig = plt.figure(figsize=(6, 6), dpi=100)
        self.three_d_bubble_ax = self.three_d_bubble_fig.add_subplot(111, projection='3d')
        self.three_d_bubble_canvas = FigureCanvasTkAgg(self.three_d_bubble_fig, self.three_d_bubble_tab_right)
        self.three_d_bubble_ax.mouse_init()

        # Jar Data
        self.load_jar_data()

        # Show Log Viewer Canvas
        self.three_d_bubble_canvas.draw()
        self.three_d_bubble_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # Event Data
        self.load_3d_bubble_data()

    def load_jar_data(self):

        # Upper Cylinder
        r = float(self.radius)
        pz = float(self.positive_z)
        # Lower Bowl/Sphere
        nz = float(self.negative_z)

        # Polar
        u = np.linspace(0, 2 * np.pi, 100)
        z = np.linspace(0, abs(pz), int((abs(pz))/2))

        U, Z = np.meshgrid(u, z)

        rstride = 1 + int((abs(pz)+abs(nz))/2/20)
        cstride = 5
        self.three_d_bubble_ax.plot_wireframe(r * cos(U), r * sin(U), np.sign(nz)*-Z, alpha=0.2, rstride=rstride, cstride=cstride)

        # Polar
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(np.pi/2, np.pi, int(abs(nz)/2))

        U, V = np.meshgrid(u, v)

        rstride = 1 + int((abs(pz)+abs(nz))/2/40)
        cstride = 5
        self.three_d_bubble_ax.plot_wireframe(r * cos(U) * sin(V), r * sin(U) * sin(V), -nz * cos(V), alpha=0.2, rstride=rstride, cstride=cstride)

        # Set Graph Limits
        f = 50
        self.three_d_bubble_ax.set_xlim(-1 * r - f, r + f)
        self.three_d_bubble_ax.set_ylim(-1 * r - f, r + f)
        if np.sign(nz) == 1:
            self.three_d_bubble_ax.set_zlim(-pz + f, 1 * nz - f)
        else:
            self.three_d_bubble_ax.set_zlim(1 * nz - f, -pz + f) 

        # Label Axes
        self.three_d_bubble_ax.set_xlabel('X (mm)')
        self.three_d_bubble_ax.set_ylabel('Y (mm)')
        self.three_d_bubble_ax.set_zlabel('Z (mm)')

    def load_3d_bubble_data(self):

        # Clear Away Previous Bubble Data
        self.three_d_bubble_ax.clear()
        self.load_jar_data()

        if self.reco_events is None:   # no reco events means we don't have reco data
            self.three_d_bubble_ax.text2D(0.05, 0.95, "No Reco Available", transform=self.three_d_bubble_ax.transAxes)
            self.three_d_bubble_canvas.draw()
            # self.logger.error('No Reco File Avaiable for Events')
            return

        if not self.show_var.get():
            self.three_d_bubble_label['state'] = tk.DISABLED
            self.three_d_bubble_position_button['state'] = tk.DISABLED
            self.show_all_reco_var.set(False)
            self.three_d_bubble_show_all_reco_checkbutton['state'] = tk.DISABLED
            self.three_d_bubble_ax.clear()
            self.load_jar_data()
            self.three_d_bubble_canvas.draw()
            return
        else:
            self.three_d_bubble_label['state'] = tk.NORMAL
            self.three_d_bubble_position_button['state'] = tk.NORMAL
            self.three_d_bubble_show_all_reco_checkbutton['state'] = tk.NORMAL

        # Set up Lists
        self.X_var = []
        self.Y_var = []
        self.Z_var = []
        self.edge_var = []

        if not self.show_all_reco_var.get():
            if not self.reco_row:  # no reco row means we don't have reco data
                self.three_d_bubble_canvas.draw()
                # self.logger.error('No Reco File Avaiable for Event')
                return
            if self.reco_row['nbub'] < 1:  # if there are no bubbles for this event
                self.three_d_bubble_canvas.draw()
                return
            for ibub in range(1, self.reco_row['nbub'] + 1):
                self.load_reco_row(ibub)

                # Load reco data
                x_var = float(self.reco_row['X'])
                if math.isnan(x_var):
                    x_var = -1000
                y_var = float(self.reco_row['Y'])
                if math.isnan(y_var):
                    y_var = -1000
                z_var = float(self.reco_row['Z'])
                if math.isnan(z_var):
                    z_var = -1000
                edge_var = float(self.reco_row['Dwall'])
                if math.isnan(edge_var):
                    edge_var = -1000

                # Assign data to variables
                self.X_var.extend([x_var])
                self.Y_var.extend([y_var])
                self.Z_var.extend([z_var])
                self.edge_var.extend([edge_var])

                # Set up Position Data Labels                
                self.bubble_label = '(%d, %d, %d)' % (x_var, y_var, z_var)
                self.edge_label = '(%.2f)' % (edge_var)

                if self.three_d_bubble_position_combobox.get() == 'Distance to Edge':
                    self.edge_label = '(%.2f)' % (edge_var)
                    # Add XYZ Position Text
                    self.three_d_bubble_ax.text(x_var, y_var, z_var, self.edge_label, alpha=0.7, color='k')

                if self.three_d_bubble_position_combobox.get() == 'XYZ Coordinates':
                    self.bubble_label = '(%d, %d, %d)' % (x_var, y_var, z_var)
                    # Add XYZ Position Text
                    self.three_d_bubble_ax.text(x_var, y_var, z_var, self.bubble_label, alpha=0.7, color='k')
                continue
            self.load_reco_row()
            
        alphavalue = 0.7
        markerSize = 15
        if self.show_all_reco_var.get():
            if self.selected_events is None:
                self.logger.error('3D Bubble Viewer: Must Apply a Cut to Show All')
                self.show_all_reco_var.set(False)
                self.load_3d_bubble_data()
                return
            else:
                self.event_length = len(self.selected_events)
                # print(self.event_length)
                if self.event_length > 1000:
                    if not tk.messagebox.askokcancel(title='Warning', message='Loading all reco bubbles may take several minutes'):
                        self.show_all_reco_var.set(False)
                        return
                    alphavalue = 0.4
                    markerSize = 10
                if self.event_length > 5000:
                    alphavalue = 0.2
                    markerSize = 5
                if self.event_length > 15000:
                    alphavalue = 0.1
                    markerSize = 4
                for selected_event in self.selected_events:
                    reco_row = selected_event
                    if not reco_row:  # no reco row means we don't have reco data
                        return
                    if reco_row['nbub'] < 1:
                        continue
                    for ibub in range(1, reco_row['nbub'] + 1):
                        # self.load_reco_row(ibub)
                        offset = ibub - 1 if ibub > 1 else 0
                        row = np.argwhere((self.reco_events['run'] == reco_row['run']) & (self.reco_events['ev'] == reco_row['ev'])).ravel()[0]
                        reco_row = self.reco_events[row + offset]

                        # Load reco data
                        x_var = float(reco_row['X'])
                        if math.isnan(x_var):
                            x_var = -1000
                        y_var = float(reco_row['Y'])
                        if math.isnan(y_var):
                            y_var = -1000
                        z_var = float(reco_row['Z'])
                        if math.isnan(z_var):
                            z_var = -1000
                        edge_var = float(reco_row['Dwall'])
                        if math.isnan(edge_var):
                            edge_var = -1000

                        # Assign data to graphed variables
                        self.X_var.extend([x_var])
                        self.Y_var.extend([y_var])
                        self.Z_var.extend([z_var])
                        self.edge_var.extend([edge_var])

                        self.bubble_label = '(%d, %d, %d)' % (x_var, y_var, z_var)
                        self.edge_label = '(%.2f)' % (edge_var)

                        if self.three_d_bubble_position_combobox.get() == 'Distance to Edge':
                            self.edge_label = '(%.2f)' % (edge_var)
                            # Add XYZ Position Text
                            self.three_d_bubble_ax.text(x_var, y_var, z_var, self.edge_label, alpha=0.7, color='k')

                        if self.three_d_bubble_position_combobox.get() == 'XYZ Coordinates':
                            self.bubble_label = '(%d, %d, %d)' % (x_var, y_var, z_var)
                            # Add XYZ Position Text
                            self.three_d_bubble_ax.text(x_var, y_var, z_var, self.bubble_label, alpha=0.7, color='k')
                    continue

        # Plot Bubble
        self.three_d_bubble_ax.scatter(self.X_var, self.Y_var, self.Z_var, alpha=alphavalue, color='r', s=markerSize)

        # Redraw Canvas with New Data
        self.three_d_bubble_canvas.draw()

    ############################################################################

    def create_three_d_bubble_widgets(self):
        
        self.three_d_bubble_tab = tk.Frame(self.notebook)
        self.notebook.add(self.three_d_bubble_tab, text='3D Bubble Map')

        #######3d Bubble Tab#######

        # First setup frames for three_d_bubble tab
        self.three_d_bubble_tab_left = tk.Frame(self.three_d_bubble_tab, bd=5, relief=tk.SUNKEN)
        self.three_d_bubble_tab_left.grid(row=0, column=0, sticky='NW')

        self.three_d_bubble_tab_right = tk.Frame(self.three_d_bubble_tab, bd=5, relief=tk.SUNKEN)
        self.three_d_bubble_tab_right.grid(row=0, column=1, sticky='NW')

        # Now with the three_d_bubble frames setup stuff
        self.three_d_bubble_show_checkbutton = tk.Checkbutton(
            self.three_d_bubble_tab_left,
            text='Turn on',
            variable=self.show_var,
            command=self.load_3d_bubble_data)
        self.three_d_bubble_show_checkbutton.grid(row=0, column=0, columnspan=2, sticky='WE')

        self.three_d_bubble_label = tk.Label(self.three_d_bubble_tab_left, text='Position Data Shown:')
        self.three_d_bubble_label.grid(row=1, column=0, sticky='WE')
        self.three_d_bubble_position_combobox = ttk.Combobox(
            self.three_d_bubble_tab_left,
            values=['None', 'XYZ Coordinates', 'Distance to Edge'])
        self.three_d_bubble_position_combobox.grid(row=1, column=1, sticky='WE')

        self.three_d_bubble_position_button = tk.Button(
            self.three_d_bubble_tab_left,
            text='Apply',
            command=self.load_3d_bubble_data)
        self.three_d_bubble_position_button.grid(row=1, column=3, sticky='NW')

        self.three_d_bubble_show_all_reco_checkbutton = tk.Checkbutton(
            self.three_d_bubble_tab_left,
            text='Show All Events for Matching Cuts in Reco File',
            variable=self.show_all_reco_var,
            command=self.load_3d_bubble_data)
        self.three_d_bubble_show_all_reco_checkbutton.grid(row=2, column=0, columnspan=2, sticky='WE')

        ############################
