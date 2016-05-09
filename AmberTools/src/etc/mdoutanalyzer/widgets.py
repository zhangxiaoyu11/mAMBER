"""
List of useful widgets we use here.
"""
from __future__ import division

from csv import writer
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
         NavigationToolbar2TkAgg)
from mdoutanalyzer.constants import LABEL_DESC
from mdoutanalyzer.windows import TextWindow
import numpy as np
from tkFileDialog import asksaveasfilename
from tkMessageBox import showerror, showwarning, showinfo
from Tkinter import *
try:
   from scipy.stats import gaussian_kde
except ImportError:
   gaussian_kde = None

def make_figure_window(fcn):
   """
   Decorator to set up a new window to draw the plot in and set up the plot
   appearance
   """
   def new_fcn(self):
      # Generate the new window
      new_window = Toplevel(self)
      fig = plt.figure(self.master.master.next_fig())
      newfig = FigureCanvasTkAgg(fig, master=new_window)
      ax = []
      # Determine our default x- and y-labels
      if self.graph_props.xlabel():
         xlab = self.graph_props.xlabel()
      else:
         xlab = self.default_xlabel()
      if self.graph_props.ylabel():
         ylab = self.graph_props.ylabel()
      else:
         ylab = self.default_ylabel()
      # Add all of the new subplots that we want
      if not self.graph_props.separate():
         nax = 1
      else:
         nax = 0
         for val in self.activelist:
            if val.get(): nax += 1
      ncols = min(nax, 3)
      nrows = nax // 3
      if nax % 3 != 0: nrows += 1
      for i in range(nax):
         newax = fig.add_subplot(nrows, ncols, i+1)
         newax.set_xlabel(xlab)
         newax.set_ylabel(ylab)
         if self.graph_props.title():
            title = self.graph_props.title(i)
         else:
            title = self.default_title(i+1)
         newax.set_title(title)
         newax.grid(self.graph_props.gridlines())
         ax.append(newax)
      # Call our decorated function and bail out if it failed
      if fcn(self, newfig, ax, self.graph_props.nexcl()) == 1:
         # This means we errored
         new_window.destroy()
         return
      # Decorate our new window and register it so it keeps a reference
      self.new_windows.append(new_window)
      newfig.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
      # Add the toolbar
      toolbar = NavigationToolbar2TkAgg(newfig, new_window)
      toolbar.update()
      newfig._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

   return new_fcn

def check_valid_input(fcn):
   """ Decorator set up to check that we have valid input """
   def new_fcn(self):
      if sum([v.get() for v in self.activelist]) == 0:
         showerror('No Data Sets!', 'No data sets chosen!', parent=self)
         return 1

      nexcl = self.graph_props.nexcl()
      if nexcl < 0 or nexcl > len(self.datasets[self.keylist[0]]):
         showerror('Bad Exclusions!',
                   'Number of excluded points must be greater than '
                   'zero and less than the number of data points!',
                   parent=self)
         return 1
      return fcn(self)

   return new_fcn

class InputEntryWindow(Toplevel):
   """ Widget that takes a single text input """
   
   def __init__(self, master, var, label, nchar=50):
      Toplevel.__init__(self, master)
      self.entry = Entry(self, textvariable=var, width=nchar)
      self.label = Label(self, text=label)
      self.button = Button(self, text='OK', command=self.destroy)
      self.entry.pack()
      self.label.pack()
      self.button.pack(fill=X)
      self.resizable(False, False)
      self.grab_set()

class LabeledEntry(Frame):
   """ This is a labeled entry widget """

   def __init__(self, master, text, *args, **kwargs):
      Frame.__init__(self, master)
      self.entry = Entry(self, *args, **kwargs)
      self.label = Label(self, text=text)
      self.entry.pack(fill=BOTH)
      self.label.pack(fill=BOTH)

class DataButton(Checkbutton):
   """ Button for each data set parsed from the mdout file """
   
   def __init__(self, master, var, label):
      Checkbutton.__init__(self, master, indicatoron=False, text=label,
                           variable=var, width=20, height=2)
      self.bind('<Button-3>', lambda event: showinfo('Energy term',
                  LABEL_DESC[label], master=self))

class _AnaButton(Button):
   """ General base class for analysis buttons """

   def __init__(self, master, datasets, activelist, graph_props, text):
      Button.__init__(self, master, text=text, width=20,
                      height=2, command=self.execute)
      self.datasets = datasets
      self.graph_props = graph_props
      self.activelist = activelist
      self.keylist = self.datasets.keys()
      self.new_windows = []

   def destroy(self):
      """ Destroy all windows we've spawned, then destroy myself """
      for window in self.new_windows:
         window.destroy()
      Button.destroy(self)

   def default_xlabel(self):
      """ Get the default label for the X-axis """
      if self.graph_props.use_time():
         return 'Time (ps)'
      else:
         return 'Frame #'

   def default_ylabel(self):
      return ""

   def default_title(self, num):
      """ Generate the default title from the data """
      nax = 0
      for v in self.activelist:
         if v.get(): nax += 1
      if not self.graph_props.separate() and nax > 1:
         return ""
      nfound = 0
      for i, key in enumerate(self.keylist):
         if self.activelist[i].get():
            nfound += 1
            if nfound == num:
               return LABEL_DESC[key]
   
class GraphButton(_AnaButton):
   """ This is the button to graph all of the selected data sets """

   def default_title(self, num):
      fpart = _AnaButton.default_title(self, num).strip()
      if fpart:
         if self.graph_props.use_time() and 'TIME(PS)' in self.keylist:
            return fpart + ' vs. Time'
         else:
            return fpart + ' vs. Frame'
      else:
         return 'Raw Data'

   @check_valid_input
   @make_figure_window
   def execute(self, fig, ax, nexcl):
      """ Graphs the data sets """
      # Try to get our x data from the Time
      try:
         xdata = self.datasets['TIME(PS)'].copy()[nexcl::
                                                  self.graph_props.stride()]
      except KeyError:
         xdata = np.arange(nexcl+1, len(self.datasets[self.keylist[0]])+1,
                           self.graph_props.stride())
      
      if not self.graph_props.use_time():
         xdata = np.arange(nexcl+1, len(self.datasets[self.keylist[0]])+1, 
                           self.graph_props.stride())
      nplt = 0
      for i, a in enumerate(self.activelist):
         if not a.get(): continue
         # plot me
         props = self.graph_props.graph_options()
         if self.graph_props.legend() and not self.graph_props.separate():
            label = self.keylist[i]
         else:
            label = '_nolegend_'
         # Catch instance where an energy is not printed on the first frame
         # (e.g., for EAMBER when restraints are on)
         ax[nplt].plot(xdata, self.datasets[self.keylist[i]].copy()[
                                          nexcl::self.graph_props.stride()],
                  label=label, **props)
         # Show the legend or not
         if self.graph_props.legend() and not self.graph_props.separate():
            ax[nplt].legend(loc=0)
         if self.graph_props.separate(): nplt += 1
      self.graph_props.reset_props()
      fig.show()

class SaveButton(_AnaButton):
   """ For saving the data to a file """

   @check_valid_input
   def execute(self):
      """ Saves the data sets to a text or CSV file """
      fname = asksaveasfilename(parent=self, defaultextension='.dat',
                                filetypes=[('Data File', '*.dat'),
                                           ('CSV File', '*.csv'),
                                           ('All Files', '*')])
      nexcl = self.graph_props.nexcl()
      xdata = np.arange(1, len(self.datasets[self.keylist[0]])+1)
      actives, keys = [xdata], ['Frame']
      for i, val in enumerate(self.activelist):
         if not val.get(): continue
         keys.append(self.keylist[i])
         actives.append(self.datasets[self.keylist[i]])
      
      f = open(fname, 'w')
      # Detect csv or not
      if fname.endswith('.csv'):
         csvwriter = writer(f)
         # Header
         csvwriter.writerow(keys)
         for i in range(nexcl, len(self.datasets[keys[1]]),
                        self.graph_props.stride()):
            csvwriter.writerow([v[i] for v in actives])
      else:
         f.write('#' + ''.join(['%16s' % n for n in keys]) + '\n')
         for i in range(nexcl, len(self.datasets[keys[1]]),
                        self.graph_props.stride()):
            f.write(' ' + ''.join(['%16.4f' % v[i] for v in actives]) + '\n')
      f.close()

class StatButton(_AnaButton):
   """ Prints average and standard deviation of each data set """

   NVALS = 4

   @check_valid_input
   def execute(self):

      nexcl = self.graph_props.nexcl()

      header = ('%-16s' + '%16s'*self.NVALS) % ('Data Set', 'Min', 'Max',
                'Avg.', 'Std. Dev.')
      report_str = header + '\n' + '-'*len(header) + '\n'

      for i, val in enumerate(self.activelist):
         if not val.get(): continue
         dset = self.datasets[self.keylist[i]][nexcl::self.graph_props.stride()]
         report_str += ('%-16s' + '%16.4f'*self.NVALS + '\n') % (self.keylist[i],
                       dset.min(), dset.max(), dset.mean(), dset.std())
      
      window = TextWindow(self.master, width=len(header)+4, height=i+5)
      window.write(report_str)
      window.grab_set()

class HistButton(_AnaButton):
   """ Histograms the data """

   def default_title(self, num):
      fpart = _AnaButton.default_title(self, num).strip()
      if fpart:
         if self.graph_props.use_kde():
            return fpart + ' Smoothed Normalized Histogram (KDE)'
         elif self.graph_props.normalize():
            return fpart + ' Normalized Histogram'
         else:
            return fpart + ' Histogram'
      elif self.graph_props.use_kde():
         return 'Smoothed Normalized Histogram (KDE)'
      elif self.graph_props.normalize():
         return 'Normalized Histograms'
      else:
         return 'Histograms'

   def default_xlabel(self):
      return ""

   def default_ylabel(self):
      if self.graph_props.normalize() or self.graph_props.kde():
         return 'Probability'
      else:
         return 'Frequency'

   @check_valid_input
   @make_figure_window
   def execute(self, fig, ax, nexcl):
      """ Graphs the histograms of the data sets """

      if self.graph_props.use_kde() and gaussian_kde is None:
         showwarning('No scipy!', 'You must have scipy in order to use a '
                     'kernel density estimate to smooth the histograms!',
                     parent=self)
         self.graph_props.noscipy()
      
      nplt = 0
      # Set the graph properties
      for i, a in enumerate(self.activelist):
         if not a.get(): continue
         dset = (self.datasets[self.keylist[i]].copy()
                        [nexcl::self.graph_props.stride()])
         props = self.graph_props.graph_options()
         if self.graph_props.legend() and not self.graph_props.separate():
            label = self.keylist[i]
         else:
            label = '_nolegend_'
         # Plot either with a KDE or not
         if self.graph_props.use_kde():
            # Use a kernel density estimate for the histogramming to provide a
            # smooth curve
            kde = gaussian_kde(dset)
            # Use 200 points to characterize the surface. Go 1/100th of the range
            # out past either side of the max and min
            kmin = dset.min() - (dset.max() - dset.min()) / 100
            kmax = dset.max() + (dset.max() - dset.min()) / 100
            xdata = np.arange(kmin, kmax+0.000000001, (kmax-kmin)/200)
            ydata = np.asarray([kde.evaluate(x) for x in xdata])
            ax[nplt].plot(xdata, ydata, label=label, **props)
         else:
            # No KDE -- straight-out histogramming
            bw = self.graph_props.binwidth()
            if bw == 0:
               bw = 3.5 * dset.std() / len(dset) ** (1/3)
            if bw > 0:
               nbins = int((np.max(dset) - np.min(dset)) / bw)
            else:
               nbins = -int(bw)
            try:
               hist, bin_edges = np.histogram(dset, nbins,
                                           density=self.graph_props.normalize())
            except TypeError:
               hist, bin_edges = np.histogram(dset, nbins,
                                           normed=self.graph_props.normalize())
            ax[nplt].plot(bin_edges[:len(hist)], hist, label=label, **props)
         # Plot our function
         if self.graph_props.legend() and not self.graph_props.separate():
            ax[nplt].legend(loc=0)
         if self.graph_props.separate(): nplt += 1
      self.graph_props.reset_props()
      fig.show()

class AutoCorrButton(_AnaButton):
   """
   Does the normed autocorrelation function (with plotting) of the data set
   """

   def default_title(self, num):
      fpart = _AnaButton.default_title(self, num).strip()
      if fpart:
         return fpart + ' Autocorrelation Function'
      else:
         return "Normalized Autocorrelation Function"

   def default_xlabel(self):
      if self.graph_props.use_time() and 'TIME(PS)' in self.keylist:
         return "Lag Time (ps)"
      return "Lag (Frame #)"
   
   def default_ylabel(self):
      return "Normalized Autocorrelation"

   @check_valid_input
   @make_figure_window
   def execute(self, fig, ax, nexcl):

      # Try to get our x data from the Time
      try:
         xdata = self.datasets['TIME(PS)'].copy()[nexcl::
                                                  self.graph_props.stride()]
      except KeyError:
         xdata = np.arange(nexcl+1, len(self.datasets[self.keylist[0]])+1,
                           self.graph_props.stride())
      
      if not self.graph_props.use_time():
         xdata = np.arange(nexcl+1, len(self.datasets[self.keylist[0]])+1, 
                           self.graph_props.stride())
      
      nplt = 0
      for i, a in enumerate(self.activelist):
         if not a.get(): continue
         # plot me
         props = self.graph_props.graph_options()
         if self.graph_props.legend() and not self.graph_props.separate():
            label = self.keylist[i]
         else:
            label = '_nolegend_'
         dset = (self.datasets[self.keylist[i]].copy()
                                 [nexcl::self.graph_props.stride()])
         dset -= dset.sum() / len(dset)
         dset /= dset.std()
         dset2 = dset.copy() / len(dset)
         acor = np.correlate(dset, dset2, 'full')
         acor = acor[len(acor)//2:]
         xend = len(acor)
         ax[nplt].plot(xdata[:xend], acor, label=label, **props)
         if self.graph_props.legend() and not self.graph_props.separate():
            ax[nplt].legend(loc=0)
         if self.graph_props.separate(): nplt += 1
      self.graph_props.reset_props()
      fig.show()

class RunningAvgButton(_AnaButton):
   """
   Plots the running average of a variable
   """

   def default_title(self, num):
      fpart = _AnaButton.default_title(self, num).strip()
      if fpart:
         return fpart + ' Running Average (window = %d)' % (
                           self.graph_props.window())
      else:
         return "Running Average (window = %d)" % (self.graph_props.window())

   @check_valid_input
   @make_figure_window
   def execute(self, fig, ax, nexcl):
      
      try:
         xdata = self.datasets['TIME(PS)'].copy()[nexcl::
                                                  self.graph_props.stride()]
      except KeyError:
         xdata = np.arange(nexcl+1, len(self.datasets[self.keylist[0]])+1,
                           self.graph_props.stride())

      if self.graph_props.use_time():
         xdata = np.arange(nexcl+1, len(self.datasets[self.keylist[0]])+1,
                           self.graph_props.stride())
      
      nplt = 0
      # Set the graph properties
      for i, a in enumerate(self.activelist):
         if not a.get(): continue
         # plot me
         props = self.graph_props.graph_options()
         if self.graph_props.legend() and not self.graph_props.separate():
            label = self.keylist[i]
         else:
            label = '_nolegend_'
         dset = (self.datasets[self.keylist[i]].copy()
                     [nexcl::self.graph_props.stride()])
         ravg = dset.copy()
         self._calc_run_avg(ravg, dset)
         ax[nplt].plot(xdata, ravg, label=label, **props)
         if self.graph_props.legend() and not self.graph_props.separate():
            ax[nplt].legend(loc=0)
         if self.graph_props.separate(): nplt += 1
      self.graph_props.reset_props()
      fig.show()
  
   def _calc_run_avg(self, ravg, dset):
      """ Calculates the running average """
      window = self.graph_props.window()
      for i in range(len(ravg)):
         ravg[i] = (np.sum(dset[max(0,i-window):min(len(ravg),i+window)]) /
                       (min(len(ravg), i+window) - max(0, i-window)))

class CumulativeAvgButton(RunningAvgButton):
   
   def default_title(self, num):
      fpart = _AnaButton.default_title(self, num).strip()
      if fpart:
         return fpart + ' Cumulative Running Average'
      else:
         return 'Cumulative Running Average'
   
   def _calc_run_avg(self, ravg, dset):
      """ Calculates the running average """
      for i in range(len(ravg)):
         ravg[i] = np.sum(dset[:i+1]) / (i+1)

if __name__ == '__main__':
   root = Tk()
   myvar = StringVar()
   ient = InputEntryWindow(root, myvar, 'Test!')
   root.mainloop()

   print 'myvar == ', myvar.get()
