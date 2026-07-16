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

        # Single mutually-exclusive view mode for the 3D tab.
        # 'off' means don't render the graph
        self.view_mode = tk.StringVar(value='off')  # 'off' | 'event' | 'run' | 'cuts'
        self._last_view_mode = None   # tracks the previously-rendered mode
        self._view_radios = {}
        self._view_base_text = {}

        self.X_var = 0
        self.Y_var = 0
        self.Z_var = 0

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

    def _event_xyz(self, reco3d, ev):
        # Representative 3D point for an event: the first non-NaN coords_3D row.
        # reco.sbc has one row per frame in frame order, so the first valid row is
        # the earliest reconstructed (bubble-formation) position. Returns None if the
        # event has no reconstruction.
        if reco3d is None:
            return None
        coords = np.asarray(reco3d['coords_3D'][reco3d['ev'] == int(ev)], dtype=float)
        if coords.size == 0:
            return None
        valid = ~np.isnan(coords).any(axis=1)
        if not valid.any():
            return None
        return coords[valid][0]

    def _add_xyz(self, xyz):
        # Append a point to the plotted lists and, if requested, label it.
        x_var, y_var, z_var = float(xyz[0]), float(xyz[1]), float(xyz[2])
        self.X_var.append(x_var)
        self.Y_var.append(y_var)
        self.Z_var.append(z_var)
        if self.three_d_bubble_position_combobox.get() == 'XYZ Coordinates':
            self.three_d_bubble_ax.text(x_var, y_var, z_var,
                                        '(%d, %d, %d)' % (x_var, y_var, z_var),
                                        alpha=0.7, color='k')

    def _refresh_view_counts(self):
        # Update the count shown next to each radio
        if not self._view_radios:
            return
        reco3d = self.reco3d_events
        ev = self.event

        # if we have valid recon, event, and a valid 3D point, the count is 1
        # FIXME: when multibubble lands we need to change this
        event_n = 1 if (reco3d is not None and ev is not None
                        and self._event_xyz(reco3d, ev) is not None) else 0
        run_n = 0
        if reco3d is not None:
            run_n = sum(self._event_xyz(reco3d, e) is not None
                        for e in np.unique(reco3d['ev']))
        if self.selected_events is None:
            cut_lbl = str(run_n)  # no cut applied means cut-matched shows the whole run
        else:
            run_events = self.selected_events[self.selected_events['run'] == self.run]
            cut_lbl = str(sum(self._event_xyz(reco3d, e['ev']) is not None
                              for e in run_events))

        self._view_radios['event'].config(text='{} ({})'.format(self._view_base_text['event'], event_n))
        self._view_radios['run'].config(text='{} ({})'.format(self._view_base_text['run'], run_n))
        self._view_radios['cuts'].config(text='{} ({})'.format(self._view_base_text['cuts'], cut_lbl))

    def load_3d_bubble_data(self):

        mode = self.view_mode.get()
        self._refresh_view_counts()

        if mode == 'off' and self._last_view_mode == 'off':
            return
        self._last_view_mode = mode

        # Clear Away Previous Bubble Data
        self.three_d_bubble_ax.clear()
        self.load_jar_data()

        if mode == 'off':
            self.three_d_bubble_label['state'] = tk.DISABLED
            self.three_d_bubble_position_button['state'] = tk.DISABLED
            self.three_d_bubble_canvas.draw()
            return
        else:
            self.three_d_bubble_label['state'] = tk.NORMAL
            self.three_d_bubble_position_button['state'] = tk.NORMAL

        # Pick which events to plot for this mode. Cut-matched is scoped to the current
        # run. With no cut applied it shows the whole run
        all_run_evs = (np.unique(self.reco3d_events['ev'])
                       if self.reco3d_events is not None else [])
        if mode == 'event':
            ev_ids = [self.event] if self.event is not None else []
        elif mode == 'cuts' and self.selected_events is not None:
            ev_ids = self.selected_events[self.selected_events['run'] == self.run]['ev']
        else:  # 'run', or 'cuts' with no active cut -> whole run
            ev_ids = all_run_evs

        # Set up Lists
        self.X_var = []
        self.Y_var = []
        self.Z_var = []
        for ev in ev_ids:
            xyz = self._event_xyz(self.reco3d_events, ev)
            if xyz is not None:
                self._add_xyz(xyz)

        # Plot Bubbles
        self.three_d_bubble_ax.scatter(self.X_var, self.Y_var, self.Z_var, alpha=0.7, color='r', s=15)

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

        tk.Label(self.three_d_bubble_tab_left, text='View:').grid(row=0, column=0, sticky='W')
        view_options = [
            ('Off', 'off'),
            ('Current event', 'event'),
            ('All events in run', 'run'),
            ('Cut-matched (this run)', 'cuts'),
        ]
        self._view_radios = {}
        self._view_base_text = {}
        for i, (text, value) in enumerate(view_options):
            rb = tk.Radiobutton(
                self.three_d_bubble_tab_left,
                text=text,
                variable=self.view_mode,
                value=value,
                command=self.load_3d_bubble_data)
            rb.grid(row=i, column=1, sticky='W')
            self._view_radios[value] = rb
            self._view_base_text[value] = text

        self.three_d_bubble_label = tk.Label(self.three_d_bubble_tab_left, text='Position Data Shown:')
        self.three_d_bubble_label.grid(row=len(view_options), column=0, sticky='WE')
        self.three_d_bubble_position_combobox = ttk.Combobox(
            self.three_d_bubble_tab_left,
            values=['None', 'XYZ Coordinates'])
        self.three_d_bubble_position_combobox.grid(row=len(view_options), column=1, sticky='WE')

        self.three_d_bubble_position_button = tk.Button(
            self.three_d_bubble_tab_left,
            text='Apply',
            command=self.load_3d_bubble_data)
        self.three_d_bubble_position_button.grid(row=len(view_options), column=3, sticky='NW')

        ############################
