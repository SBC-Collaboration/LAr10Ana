# Imports
import os
import re
import matplotlib
import numpy as np
from pylab import *
import pandas as pd
import tkinter as tk
from glob import glob
#
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from scipy.optimize import curve_fit
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import sys


class LogViewer(tk.Frame):
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)

        # Log Viewer Variables
        self.fitlines = []
        self.datalines = []
        self.y_headers_right = []
        self.log_fits = ['Linear', 'Quadratic', 'Gaussian', 'Exponential', 'Logarithmic']

        # Initial Functions
        self.populate_logs()
        self.create_logviewer_widgets()
        self.logviewer_canvas_setup()

        self.xvarsel_prev = self.xvarsel.get()
        self.prev_chosen_log = self.logvarsel.get()

        # Hide Log Viewer Widgets
        self.bottom_frame_5.grid_remove()
        self.fra2.pack_forget()

    def logviewer_canvas_setup(self):
        # These used to be created every time a new event was loaded or something
        # changed in the log viewer. This lead to many figs, ax, and canvas that
        # continued to fill memory until python crashed.

        # Create Figures, Axes, and Canvases
        self.log_viewer_fig = Figure(figsize=(13, 5), dpi=100)
        self.fig_width = self.log_viewer_fig.get_size_inches()[0] * self.log_viewer_fig.dpi
        self.log_viewer_ax = self.log_viewer_fig.add_subplot(111)
        self.log_viewer_canvas = FigureCanvasTkAgg(self.log_viewer_fig, self.graph)

        # Show Log Viewer Canvas
        self.log_viewer_canvas.draw()
        self.log_viewer_canvas.get_tk_widget().grid(row=0, column=1)

        # Create Log Viewer Toolbar
        self.nav = NavigationToolbar2Tk(self.log_viewer_canvas, self.bottom_frame_5)
        self.nav.pack(side='top')

    #################################LOG VIEWER#################################

    def populate_logs(self):
        # Get List of Logs
        path = os.path.join(self.log_directory, '2*')
        self.list_logs = glob(path)

        # Logs Not Found
        if not self.list_logs:
            #             logger.error('Could not find log files. Check that the directory is correct:\n{}'.format(path))
            self.logger.warning('Could not find log files. Check that the directory is correct:\n{}'.format(path))
            self.list_logs = ' '
            return

        # Clean Up String For Visibility
        convention = re.compile('\d+_\d+')
        self.list_logs = [convention.findall(folder) for folder in self.list_logs]

        # Reformat to Single List of Strings
        for i in range(0, len(self.list_logs)):
            self.list_logs[i] = self.list_logs[i][0]

    def load_log_data(self, error_reset):    
        print('in logviewer, logpath: ', self.log_directory)
        # Get Currently Selected Log
        chosen_log = self.logvarsel.get()

        # Stop If Invalid Log File or Same Log File
        if (
                chosen_log == ' ' or chosen_log == 'Choose a Log File' or chosen_log == self.prev_chosen_log) and not error_reset:
            self.prev_chosen_log = chosen_log
            return
        else:
            self.prev_chosen_log = chosen_log

        # When an error occurs, no need to load file, data, x_variable menu again
        if not error_reset:
            # File Path
            file = os.path.join(self.log_directory, chosen_log, '{}.txt'.format(chosen_log))

            # Read Headers From File
            with open(file, 'r') as f:
                row_with_headers = 1
                for row_num, row_data in enumerate(f):
                    if row_num == row_with_headers:
                        self.x_headers = row_data.split()
                        break

            # Load Data
            first_row_with_data = 7
            self.data = pd.read_csv(file, skiprows=first_row_with_data - 1, header=None, delim_whitespace=True)

            # Don't Reset If Already Chosen
            if self.xvarsel.get() == 'Choose an X-Variable' or self.xvarsel.get() == '':
                # Destroy Old Dropdown
                self.x_select.destroy()

                # Select X Variable
                self.xvarsel = tk.StringVar()
                self.xvarsel.set('Choose an X-Variable')
                self.x_select = tk.OptionMenu(self.fra3, self.xvarsel, *self.x_headers, command=self.update_data)
                self.x_select.pack(side='top')

        # Load Y Headers
        self.y_headers_left = self.x_headers[:]

        # Remove Session_ID From List - Not Needed and Causes Errors
        try:
            index = self.y_headers_left.index('sessionID')
            del self.y_headers_left[index]
        except:
            pass

        try:
            # Remove Any Element Already Plotted
            for item in self.y_headers_right:
                # Find Index in Tracker Array
                rowindex = [row[0] for row in self.datalines].index(item)

                # Remove If Same Log File
                if self.logvarsel.get() == self.datalines[rowindex][2]:
                    # Split to Y-Variable Name
                    itemsplit = item.split(' ', 1)[0]

                    # Remove From Left Side
                    self.y_headers_left.remove(itemsplit)
        except Exception as e:
            print(e, 'Line: ', sys.exc_info()[2].tb_lineno)

        # Clear Old Y Listbox
        self.listbox_not_sel.delete(0, tk.END)

        # Populate Y Listbox
        for element in self.y_headers_left:
            self.listbox_not_sel.insert(tk.END, element)

    def temp_log_data(self, log):
        # File Path
        file = os.path.join(self.log_directory, log, '{}.txt'.format(log))

        # Read Headers From File
        with open(file, 'r') as f:
            row_with_headers = 1
            for row_num, row_data in enumerate(f):
                if row_num == row_with_headers:
                    new_x_headers = row_data.split()
                    break

        # Load Data
        first_row_with_data = 7
        new_data = pd.read_csv(file, skiprows=first_row_with_data - 1, header=None, delim_whitespace=True)

        return new_x_headers, new_data

    def show_plot(self, empty):
        # Legend
        if empty:
            self.log_viewer_ax.legend_.remove()
        else:
            self.leg = self.log_viewer_ax.legend(loc='upper left', bbox_to_anchor=(1, 1))

        # Lock to Tight View
        self.log_viewer_fig.set_tight_layout(True)

        # Update Canvas
        self.log_viewer_canvas.draw()
        self.log_viewer_canvas.get_tk_widget().grid(row=0, column=1)

        # Unlock Layout
        self.log_viewer_fig.set_tight_layout(False)

        # Update Canvas to Fit Legend
        if not empty:
            # Width of Legend
            leg_width = self.leg.get_frame().get_width()

            # Location of Legend as Fraction of Total Width
            buffer = 0.01
            leg_loc = ((self.fig_width - leg_width) / self.fig_width) - buffer

            # Adjust Right Side to Fit Legend
            self.log_viewer_fig.subplots_adjust(right=leg_loc)

        # Rescale Axes
        self.log_viewer_ax.relim()
        self.log_viewer_ax.autoscale_view()

        # Update Canvas
        self.log_viewer_canvas.draw()
        self.log_viewer_canvas.get_tk_widget().grid(row=0, column=1)

    ##Fit Functions##
    def gaussian(self, x, a, x0, sigma):
        return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

    def ln(self, x, a, b):
        return a * np.log(x) + b

    def ex(self, x, a, b):
        return a * np.exp(b * x)

    #################

    def add_fitline(self, fit_type, xdata, ydata, yvarunsplit, new, lines):
        try:
            # Fit Line
            if fit_type == 'Linear':
                # Get Fit Coefficients
                fit = np.polyfit(xdata, ydata, deg=1)

                # Calculate New Y-Data
                newydata = fit[0] * xdata + fit[1]

                # Plot
                if new:
                    self.log_viewer_ax.plot(xdata, newydata,
                                            label='Linear Fit {0:s}\n {1:.4g}x + {2:.4g}'.format(yvarunsplit, fit[0],
                                                                                                 fit[1]))
                else:
                    # Set New Arrays to Plot
                    lines.set_data(xdata, newydata)

                    # Update Legend
                    lines.set_label('Linear Fit {0:s}\n {1:.4g}x + {2:.4g}'.format(yvarunsplit, fit[0], fit[1]))
            elif fit_type == 'Quadratic':
                # Get Fit Coefficients
                fit = np.polyfit(xdata, ydata, deg=2)

                # Calculate New Y-Data
                newydata = fit[0] * xdata ** 2 + fit[1] * xdata + fit[2]

                # Plot
                if new:
                    self.log_viewer_ax.plot(xdata, newydata,
                                            label='Quadratic Fit {0:s}\n {1:.4g}x$^2$ + {2:.4g}x + {3:.4g}'.format(
                                                yvarunsplit, fit[0], fit[1], fit[2]))
                else:
                    # Set New Arrays to Plot
                    lines.set_data(xdata, newydata)

                    # Update Legend
                    lines.set_label(
                        'Quadratic Fit {0:s}\n {1:.4g}x$^2$ + {2:.4g}x + {3:.4g}'.format(yvarunsplit, fit[0], fit[1],
                                                                                         fit[2]))
            elif fit_type == 'Gaussian':
                # Calculate Statistical Variables Needed
                mean = sum(xdata * ydata) / sum(ydata)
                sigma = np.sqrt(sum(ydata * (xdata - mean) ** 2) / sum(ydata))

                # Get Fit Coefficients
                popt, pcov = curve_fit(self.gaussian, xdata, ydata, p0=[max(ydata), mean, sigma])

                # Plot
                if new:
                    self.log_viewer_ax.plot(xdata, self.gaussian(xdata, *popt),
                                            label='Gaussian Fit {}'.format(yvarunsplit))
                else:
                    # Set New Arrays to Plot
                    lines.set_data(xdata, self.gaussian(xdata, *popt))

                    # Update Legend
                    lines.set_label('Gaussian Fit {}'.format(yvarunsplit))
            elif fit_type == 'Exponential':
                # Get Fit Coefficients
                popt, pcov = curve_fit(self.ex, xdata, ydata)

                # Calculate New Y-Data
                newydata = self.ex(xdata, *popt)

                # Plot
                if new:
                    self.log_viewer_ax.plot(xdata, newydata,
                                            label='Exponential Fit {0:s}\n ({1:.4g})e$^{{({2:.4g})x}}$'.format(
                                                yvarunsplit, popt[0], popt[1]))
                else:
                    # Set New Arrays to Plot
                    lines.set_data(xdata, newydata)

                    # Update Legend
                    lines.set_label(
                        'Exponential Fit {0:s}\n ({1:.4g})e$^{{({2:.4g})x}}$'.format(yvarunsplit, popt[0], popt[1]))
            elif fit_type == 'Logarithmic':
                # Store X-Data to be Modified For Calculations
                xdata_calc = xdata

                # Initial Shift
                c = 0

                # Shift Data if Negative
                try:
                    # Shift By Most Negative Number
                    c = abs(min([n for n in xdata if n < 0]))
                    xdata_calc = [n + c for n in xdata]
                except:
                    pass

                # Convert Back to Numpy Array
                xdata_calc = np.asarray(xdata_calc, dtype=float)

                # Set Indices With Zero Values to Small Numbers
                xdata_calc[xdata_calc == 0] = 1E-10

                # Get Fit Coefficients
                popt, pcov = curve_fit(self.ln, xdata_calc, ydata)

                # Calculate New Y-Data
                newydata = self.ln(xdata_calc, *popt)

                # Plot
                if new:
                    self.log_viewer_ax.plot(xdata, newydata,
                                            label='Logarithmic Fit {0:s}\n ({1:.4g})ln(x + {2:.4g}) + {3:.4g}'.format(
                                                yvarunsplit, popt[0], c, popt[1]))
                else:
                    # Set New Arrays to Plot
                    lines.set_data(xdata, newydata)

                    # Update Legend
                    lines.set_label(
                        'Logarithmic Fit {0:s}\n ({1:.4g})ln(x + {2:.4g}) + {3:.4g}'.format(yvarunsplit, popt[0], c,
                                                                                            popt[1]))
            else:
                raise Exception('Select a fit type.')
        except:
            raise Exception('An error occurred while fitting data for {}.'.format(yvarunsplit))

    def update_data(self, _):
        # Stop If Nothing Is Plotted or Same X-Variable
        if len(self.datalines) == 0 or self.xvarsel.get() == self.xvarsel_prev:
            self.xvarsel_prev = self.xvarsel.get()
            return
        else:
            self.xvarsel_prev = self.xvarsel.get()

        # Update Plotted Data and Fit Lines
        try:
            # Current Selected Log
            chosen_log = self.logvarsel.get()

            # Get Current Lines Plotted
            cur_lines_plot = self.log_viewer_ax.lines[:]

            # Update Lines to New X-Data
            for lines in cur_lines_plot:
                # Check If Line is Fitline or Dataline
                if lines in (row[2] for row in self.fitlines):
                    # Find Index in Tracker Array
                    row_index = [row[2] for row in self.fitlines].index(lines)

                    # Load Data
                    if chosen_log != self.fitlines[row_index][6]:
                        x_headers, data = self.temp_log_data(self.fitlines[row_index][6])
                    else:
                        data = self.data
                        x_headers = self.x_headers

                    # New Selected X-Variable
                    xvar = self.xvarsel.get()

                    # Relabel X-Axis
                    self.log_viewer_ax.set_xlabel(xvar)

                    # New Data
                    xindex = x_headers.index(xvar)
                    xdata = data.iloc[:, xindex].values

                    # Get Y-Variable
                    yvarunsplit = self.fitlines[row_index][0]
                    yvar = yvarunsplit.split(' ', 1)[0]

                    # Get Location of Header in Original File
                    yindex = x_headers.index(yvar)

                    # Store Column Data
                    ydata = data.iloc[:, yindex]

                    # Get Indices
                    first_index = self.fitlines[row_index][3]
                    last_index = self.fitlines[row_index][4]

                    # Trim X-Data
                    fitxdata = xdata[first_index:last_index]

                    # Trim Y-Data
                    ydata = ydata[first_index:last_index]

                    # Get Type of Fit
                    fit_type = self.fitlines[row_index][5]

                    # Prevent Improper Fitting
                    if len(set(xdata)) == 1:
                        print('Infinite slope for {}, and cannot be fit to.'.format(yvarunsplit))
                        self.delete_fit(True, yvarunsplit)
                    else:
                        # Update Fit Line
                        try:
                            self.add_fitline(fit_type, fitxdata, ydata, yvarunsplit, False, lines)
                        except Exception as e:
                            print(e, 'Line: ', sys.exc_info()[2].tb_lineno)
                            self.delete_fit(True, yvarunsplit)
                else:
                    # Find Index in Tracker Array
                    row_index = [row[3] for row in self.datalines].index(lines)

                    # Load Data
                    if chosen_log != self.datalines[row_index][2]:
                        x_headers, data = self.temp_log_data(self.datalines[row_index][2])
                    else:
                        data = self.data
                        x_headers = self.x_headers

                    # New Selected X-Variable
                    xvar = self.xvarsel.get()

                    # Relabel X-Axis
                    self.log_viewer_ax.set_xlabel(xvar)

                    # New Data
                    xindex = x_headers.index(xvar)
                    xdata = data.iloc[:, xindex].values

                    # Update X-Data
                    lines.set_xdata(xdata)

            # Update Plot
            self.show_plot(False)
        except Exception as e:
            # Print Error
            print(e, 'Line: ', sys.exc_info()[2].tb_lineno)

            # Reset Variables
            self.fitlines = []
            self.datalines = []
            self.y_headers_right = []
            self.load_log_data(True)

            # Clear Plot
            self.log_viewer_ax.lines = []

            # Clear Right Listbox
            self.listbox_sel.delete(0, tk.END)

            # Update Plot
            self.show_plot(True)

    def populate_plot(self, y_header):
        # Selected Variables
        xvar = self.xvarsel.get()

        # Label X-Axis
        self.log_viewer_ax.set_xlabel(xvar)

        # Get X Data
        xindex = self.x_headers.index(xvar)
        xdata = self.data.iloc[:, xindex].values

        # Get Location of Header in Original File
        yindex = self.x_headers.index(y_header)

        # Store Column Data
        ydata = self.data.iloc[:, yindex]

        # Plot
        self.log_viewer_ax.plot(xdata, ydata, label='{} - {}'.format(y_header, self.logvarsel.get()))

        # Keep Track of Which Data Lines Goes With Which Variable
        self.datalines.append(
            ['{} - {}'.format(y_header, self.logvarsel.get()), len(self.log_viewer_ax.lines) - 1, self.logvarsel.get(),
             self.log_viewer_ax.lines[len(self.log_viewer_ax.lines) - 1]])

        # Update Plot
        self.show_plot(False)

    def create_trendline(self):
        # Select Y-Variable Before Fitting
        if self.yvarfit.get() == 'Choose a Y-Variable to Fit':
            print('Select a y-variable to fit first.')
            return
        elif self.yvarfit.get() == '':
            print('Cannot fit to empty dataset.')
            return

        # Selected X Variable
        xvar = self.xvarsel.get()

        # Selected Y Variable
        yvarunsplit = self.yvarfit.get()
        yvar = yvarunsplit.split(' ', 1)[0]
        yvarrun = yvarunsplit.split(' ', 2)[2]

        # Currently Selected Log
        chosen_log = self.logvarsel.get()

        # Get Data From Appropriate Log File
        if chosen_log != yvarrun:
            x_headers, data = self.temp_log_data(yvarrun)
        else:
            x_headers = self.x_headers
            data = self.data

        # Get Indices
        if self.text_box_first.get() == '' or self.text_box_first.get() == '0':
            first_index = 0
        else:
            try:
                first_index = int(self.text_box_first.get())
            except:
                print('First index is invalid (not an integer or blank).\n')
                return
        if self.text_box_last.get() == '' or self.text_box_last.get() == '0':
            last_index = len(data.iloc[:, 0])
        else:
            try:
                last_index = int(self.text_box_last.get()) + 1
            except:
                print('Last index is invalid (not an integer or blank).\n')
                return

        # Catch More Errors
        if first_index >= last_index - 1:
            print('Last index must be strictly greater than the first index.\n')
            return
        elif first_index < 0 or last_index < 0:
            print('Indices must be a positive integer or blank.\n')
            return
        elif first_index > len(data.iloc[:, 0]) or last_index > len(data.iloc[:, 0]):
            print('Index out of data range.\n')
            return

        # Get X Data
        xindex = x_headers.index(xvar)
        xdata = data.iloc[:, xindex].values

        # Trim X Data
        xdata = xdata[first_index:last_index]

        # Prevent Improper Fitting
        if len(set(xdata)) == 1:
            print('Infinite slope for {}, and cannot be fit to.'.format(yvarunsplit))
            return

        # If Fitting Variable Already Fit, Delete Plot
        if yvarunsplit in [row[0] for row in self.fitlines]:
            self.delete_fit(False, [])

        # Get Location of Header in Original File
        yindex = x_headers.index(yvar)

        # Store Column Data
        ydata = data.iloc[:, yindex].values

        # Trim Y Data
        ydata = ydata[first_index:last_index]

        # Get Type of Fit
        fit_type = self.trendlines.get()

        # Add Fit
        try:
            self.add_fitline(fit_type, xdata, ydata, yvarunsplit, True, '')
        except Exception as e:
            print(e, 'Line: ', sys.exc_info()[2].tb_lineno)
            return

        # Keep Track of Which Fit Lines Goes With Which Variable
        self.fitlines.append(['{} - {}'.format(yvar, yvarrun), len(self.log_viewer_ax.lines) - 1,
                              self.log_viewer_ax.lines[len(self.log_viewer_ax.lines) - 1], first_index, last_index,
                              fit_type, yvarrun])

        # Update Plot
        self.show_plot(False)

    def delete_fit(self, update_elsewhere, item):
        # Select Y-Variable Before Deleting
        if self.yvarfit.get() == 'Choose a Y-Variable to Fit' and not update_elsewhere:
            print('Select a y-variable to delete first.')
            return

        # Get Fit Variable
        if not item:
            yvar = self.yvarfit.get()
        else:
            yvar = item

        # Plot Fitline Before Deleting
        if yvar not in [row[0] for row in self.fitlines] and not update_elsewhere:
            print('Fitline is not currently plotted.')
            return

        # Get Index Where Variable Fit is Stored
        rowind = [row[0] for row in self.fitlines].index(yvar)

        # Get Index of Line in ax.lines
        fitindex = self.fitlines[rowind][1]

        # No Longer Plotted, Remove From Tracker
        del self.fitlines[rowind]

        # Delete Line
        del self.log_viewer_ax.lines[fitindex]

        # Shift Stored Indices Down 1
        for row in self.datalines:
            if row[1] > fitindex:
                row[1] -= 1
        for row in self.fitlines:
            if row[1] > fitindex:
                row[1] -= 1

        # Update Plot If Not Updating With Elsewhere
        if not update_elsewhere:
            self.show_plot(False)

    def delete_data(self, item):
        # Get Fit Variable
        yvar = item

        # Get Index Where Variable Fit is Stored
        rowind = [row[0] for row in self.datalines].index(yvar)

        # Get Index of Line in ax.lines
        fitindex = self.datalines[rowind][1]

        # No Longer Plotted, Remove From Tracker
        del self.datalines[rowind]

        # Delete Line
        del self.log_viewer_ax.lines[fitindex]

        # Shift Stored Indices Down 1
        for row in self.datalines:
            if row[1] > fitindex:
                row[1] -= 1
        for row in self.fitlines:
            if row[1] > fitindex:
                row[1] -= 1

        # Update Plot
        if len(self.datalines) == 0:
            self.show_plot(True)
        else:
            self.show_plot(False)

    def move_right(self):
        # Prevent From Running With Empty Selection
        if not self.listbox_not_sel.curselection():
            return

        # Prevent Moving If No X-Value
        if self.xvarsel.get() == 'Choose an X-Variable':
            print('Choose an x-variable before plotting.')
            return

        # Get Selected Item
        item = self.listbox_not_sel.get(self.listbox_not_sel.curselection())

        # Add to Graph
        self.populate_plot(item)

        # Delete From Left Side
        index = self.y_headers_left.index(item)
        del self.y_headers_left[index]
        self.listbox_not_sel.delete(index)

        # Add to Right Side
        self.y_headers_right.append('{} - {}'.format(item, self.logvarsel.get()))
        self.listbox_sel.insert(tk.END, '{} - {}'.format(item, self.logvarsel.get()))

        # Destroy Old OptionMenus
        self.yvarfit_box.destroy()

        # Update Fit Dropdown
        self.yvarfit = tk.StringVar()
        self.yvarfit.set('Choose a Y-Variable to Fit')
        self.yvarfit_box = tk.OptionMenu(self.fra8, self.yvarfit, *self.y_headers_right)
        self.yvarfit_box.pack(side='top')

    def move_left(self):
        # Prevent From Running With Empty Selection
        if not self.listbox_sel.curselection():
            return

        # Get Selected Item
        item = self.listbox_sel.get(self.listbox_sel.curselection())

        # Find Index in Tracker Array
        rowindex = [row[0] for row in self.datalines].index(item)

        # Remove Fit Line if Exists
        try:
            self.delete_fit(True, item)
        except:
            pass

        # Delete From Right Side
        index = self.y_headers_right.index(item)
        del self.y_headers_right[index]
        self.listbox_sel.delete(index)

        # Add to Left Side - Not an error reset, but want same code to run in load_log_data
        if self.logvarsel.get() == self.datalines[rowindex][2]:
            self.load_log_data(True)

        # Remove From Graph
        self.delete_data(item)

        # Destroy Old OptionMenus
        self.yvarfit_box.destroy()

        # Update Fit Dropdown - Glitches if y_headers_right is empty and use asterisk to fill
        self.yvarfit = tk.StringVar()
        self.yvarfit.set('Choose a Y-Variable to Fit')
        if self.y_headers_right:
            self.yvarfit_box = tk.OptionMenu(self.fra8, self.yvarfit, *self.y_headers_right)
        else:
            self.yvarfit_box = tk.OptionMenu(self.fra8, self.yvarfit, [])
        self.yvarfit_box.pack(side='top')

    ############################################################################

    def create_logviewer_widgets(self):
        self.log_tab = tk.Frame(self.notebook)
        self.notebook.add(self.log_tab, text='Log Viewer')

        #######Log Viewer Tab#######
        # Graph
        self.graph = tk.Frame(self.log_tab, bd=5, relief=tk.SUNKEN)
        self.graph.pack(side='top')

        # Bottom Frame - All Buttons, Boxes, Etc.
        self.fra2 = tk.Frame(self.log_tab)
        self.fra2.pack(side='top')

        # Subframe For Listboxes, Pulldown Menus
        self.fra3 = tk.Frame(self.fra2, bd=5, relief=tk.SUNKEN)
        self.fra3.pack(fill='y', side='left')

        # Subframe For Listboxes and Listbox Text
        self.fra4 = tk.Frame(self.fra2, bd=5, relief=tk.SUNKEN)
        self.fra4.pack(fill='y', side='left')

        # Subframe For Only Listbox Text
        self.fra5 = tk.Frame(self.fra4)
        self.fra5.pack(side='top')

        # Subframe For Only Listboxes
        self.fra6 = tk.Frame(self.fra4)
        self.fra6.pack(side='top')

        # Subframe For Applying Trends
        self.fra7 = tk.Frame(self.fra2, bd=5, relief=tk.SUNKEN)
        self.fra7.pack(fill='y', side='left')

        # Subframe for Selecting Y Variable to Fit
        self.fra8 = tk.Frame(self.fra7)
        self.fra8.pack(side='top')

        # Subframe For First Index
        self.fra9 = tk.Frame(self.fra7)
        self.fra9.pack(side='top')

        # Subframe For End Index
        self.fra10 = tk.Frame(self.fra7)
        self.fra10.pack(side='top')

        # Subframe For Trendline Dropdown
        self.fra11 = tk.Frame(self.fra7)
        self.fra11.pack(side='top')

        # Subframe For Delete/Apply Trendline
        self.fra12 = tk.Frame(self.fra7)
        self.fra12.pack(side='top')

        # Subframe For Horizontal Scrollbars
        self.fra13 = tk.Frame(self.fra4)
        self.fra13.pack(fill='x', side='top')

        # Y-Fit Label
        self.y_fit_label = tk.Label(self.fra8, text='Select Y-Variable to Fit:')
        self.y_fit_label.pack(side='top')

        # Selecting Y Variable to Fit Dropdown
        self.yvarfit = tk.StringVar()
        self.yvarfit.set('Choose a Y-Variable to Fit')
        self.yvarfit_box = tk.OptionMenu(self.fra8, self.yvarfit, [])
        self.yvarfit_box.pack(side='top')

        # Trim Label
        self.Trim = tk.Label(self.fra9, text='Trim Data:')
        self.Trim.pack(side='top')

        # First Index Text
        self.list_label = tk.Label(self.fra9, text='Select First Index:')
        self.list_label.pack(side='left')

        # First Index Box
        self.text_box_first = tk.Entry(self.fra9)
        self.text_box_first.pack(side='left')

        # End Index Text
        self.list_label = tk.Label(self.fra10, text='Select End Index:')
        self.list_label.pack(side='left')

        # End Index Box
        self.text_box_last = tk.Entry(self.fra10)
        self.text_box_last.pack(side='left')

        # Fit Type Label
        self.fit_type_label = tk.Label(self.fra11, text='Select Fit Type:')
        self.fit_type_label.pack(side='top')

        # Trendline Dropdown
        self.trendlines = tk.StringVar()
        self.trendlines.set('Choose a Fit Type')
        self.trendlines_sel = tk.OptionMenu(self.fra11, self.trendlines, *self.log_fits)
        self.trendlines_sel.pack(side='top')

        # Delete Trendline
        self.delete_trendline = tk.Button(self.fra12, text='Delete Trendline',
                                          command=lambda: self.delete_fit(False, []))
        self.delete_trendline.pack(side='left')

        # Apply Trendline
        self.apply_trendline = tk.Button(self.fra12, text='Apply Trendline', command=self.create_trendline)
        self.apply_trendline.pack(side='left')

        # Log Files Label
        self.log_label = tk.Label(self.fra3, text='Select Log File:')
        self.log_label.pack(side='top')

        # Select Log Files
        self.logvarsel = tk.StringVar()
        self.logvarsel.set('Choose a Log File')
        self.log_select = tk.OptionMenu(self.fra3, self.logvarsel, *self.list_logs,
                                        command=lambda _: self.load_log_data(False))
        self.log_select.pack(side='top')

        # X-Variables Label
        self.x_label = tk.Label(self.fra3, text='Select X-Variable:')
        self.x_label.pack(side='top')

        # Select X-Variables
        self.xvarsel = tk.StringVar()
        self.xvarsel.set('Choose an X-Variable')
        self.x_select = tk.OptionMenu(self.fra3, self.xvarsel, '', command=self.update_data)
        self.x_select.pack(side='top')

        # Listbox Labels
        self.list_label1 = tk.Label(self.fra5, text='Select Y-Variables:')
        self.list_label1.pack(side='left')
        self.list_label2 = tk.Label(self.fra5, text='              ')
        self.list_label2.pack(side='left')
        self.list_label3 = tk.Label(self.fra5, text='Y-Variables Plotted:')
        self.list_label3.pack(side='left')

        # Listbox Y-Variable - Not Selected
        self.scrollbar_left_vert = tk.Scrollbar(self.fra6, orient='vertical')
        self.scrollbar_left_hori = tk.Scrollbar(self.fra13, orient='horizontal')
        self.listbox_not_sel = tk.Listbox(self.fra6, xscrollcommand=self.scrollbar_left_hori.set,
                                          yscrollcommand=self.scrollbar_left_vert.set)
        self.scrollbar_left_vert.config(command=self.listbox_not_sel.yview)
        self.scrollbar_left_hori.config(command=self.listbox_not_sel.xview)
        self.listbox_not_sel.pack(side='left')
        self.scrollbar_left_vert.pack(side='left', fill='y')
        self.scrollbar_left_hori.grid(row=0, column=0, sticky=tk.N + tk.S + tk.W + tk.E)
        self.listbox_not_sel.bind('<<ListboxSelect>>', lambda _: self.move_right())

        # Listbox Y-Variable - Selected
        self.scrollbar_right_vert = tk.Scrollbar(self.fra6, orient='vertical')
        self.scrollbar_right_hori = tk.Scrollbar(self.fra13, orient='horizontal')
        self.listbox_sel = tk.Listbox(self.fra6, xscrollcommand=self.scrollbar_right_hori.set,
                                      yscrollcommand=self.scrollbar_right_vert.set)
        self.scrollbar_right_vert.config(command=self.listbox_sel.yview)
        self.scrollbar_right_hori.config(command=self.listbox_sel.xview)
        self.listbox_sel.pack(side='left')
        self.scrollbar_right_vert.pack(side='left', fill='y')
        self.scrollbar_right_hori.grid(row=0, column=1, sticky=tk.N + tk.S + tk.W + tk.E)
        self.listbox_sel.bind('<<ListboxSelect>>', lambda _: self.move_left())

        # Horizontal Scrollbars Fill Frame Evenly
        for i in range(2):
            tk.Grid.columnconfigure(self.fra13, i, weight=1)

        ############################
