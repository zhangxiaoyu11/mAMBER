"""
Data structure controlling the appearance of the graph and a GUI control panel
exposing those options
"""

from mdoutanalyzer.widgets import LabeledEntry
from Tkinter import *

class GraphProperties(object):
   """ Defines the graph properties """

   def __init__(self):
      self._xlabel = StringVar()
      self._ylabel = StringVar()
      self._title = StringVar()
      self._linewidth = DoubleVar()
      self._linestyle = StringVar()
      self._pointstyle = StringVar()
      self._pointsize = DoubleVar()
      self._legend = BooleanVar()
      self._bin_width = DoubleVar()
      self._normalize = BooleanVar()
      self._use_kde = BooleanVar()
      self._gridlines = BooleanVar()
      self._nexcl = IntVar()
      self._stride = IntVar()
      self._window = IntVar()
      self._use_time = BooleanVar()
      self._separate = BooleanVar()
      # Set defaults
      self._legend.set(True)
      self._normalize.set(True)
      self._use_kde.set(False)
      self._gridlines.set(True)
      self._use_time.set(True)
      self._separate.set(False)
      self._bin_width.set(0.0)
      self._nexcl.set(0)
      self._stride.set(1)
      self._window.set(10)
      self._linewidth.set(1.0)
      self._pointsize.set(6.0)
      self._linestyle.set('-')
      self._pointstyle.set('')
      self._colorschemes = ['k', 'b', 'r', 'g', 'c', 'm', 'y']
      self._pointstyles = ['.', ',', 'o', 'v', '^', '<', '>', 's', '*',
                           'p', 'h', 'H', '+', 'D']
      self._linestyles = ['-', '--', '-.', ':']
      self._used_lines = 0
      self._used_points = 0
      self._used_colors = 0
   
   def lines(self):
      return self.linestyle().strip() == ''

   def points(self):
      return self.pointstyle().strip() == ''

   def xlabel(self):
      return self._xlabel.get()

   def ylabel(self):
      return self._ylabel.get()

   def title(self):
      return self._title.get()

   def legend(self):
      return self._legend.get()

   def get_next_color(self):
      i = self._used_colors
      self._used_colors += 1
      return self._colorschemes[i % len(self._colorschemes)]

   def reset_props(self):
      self._used_colors = 0
      self._used_points = 0
      self._used_lines = 0

   def linewidth(self):
      return self._linewidth.get()

   def binwidth(self):
      return self._bin_width.get()

   def normalize(self):
      return self._normalize.get()

   def gridlines(self):
      return bool(self._gridlines.get())

   def use_time(self):
      return self._use_time.get()

   def linestyle(self):
      return self._linestyle.get()

   def use_kde(self):
      return self._use_kde.get()

   def pointstyle(self):
      return self._pointstyle.get()
   
   def pointsize(self):
      return self._pointsize.get()

   def nexcl(self):
      return self._nexcl.get()

   def stride(self):
      return self._stride.get()

   def window(self):
      return self._window.get()

   def separate(self):
      return self._separate.get()

   def graph_options(self):
      """ Returns a dict of the graph options """
      r = {'linestyle' : self.linestyle(),
           'marker' : self.pointstyle(),
           'linewidth' : self.linewidth(),
           'color' : self.get_next_color(),
           'markersize' : self.pointsize()
          }
      if self.linestyle() == 'u':
         r['linestyle'] = self._linestyles[self._used_lines %
                                           len(self._linestyles)]
         self._used_lines += 1
      if self.pointstyle() == 'u':
         r['marker'] = self._pointstyles[self._used_points %
                                         len(self._pointstyles)]
         self._used_points += 1
      
      return r
   
   def noscipy(self):
      """ Call this if we have no scipy -- turns off scipy-only options """
      self._use_kde.set(False)


class GraphControlWindow(Toplevel):
   """ Frame controlling the graph options """

   def __init__(self, master, graph_props):
      Toplevel.__init__(self, master)
      self.title('Graph Property Control')
      self.check_frame = Frame(self)
      self.resizable(False, False)
      self.graph_props = graph_props
      # Draw the checkbuttons
      self.legend = Checkbutton(self.check_frame, text='Show Legend',
                       onvalue=True, offvalue=False,
                       variable=self.graph_props._legend)
      self.gridl = Checkbutton(self.check_frame, text='Show Grid Lines',
                       onvalue=True, offvalue=False,
                       variable=self.graph_props._gridlines)
      self.norm = Checkbutton(self.check_frame, text='Normalize Histograms',
                       onvalue=True, offvalue=False,
                       variable=self.graph_props._normalize)
      self.time = Checkbutton(self.check_frame, text='Use Time as X-Coord',
                       onvalue=True, offvalue=False,
                       variable=self.graph_props._use_time)
      self.kde = Checkbutton(self.check_frame,
                       text='Smooth Histograms (KDE)', onvalue=True,
                       offvalue=False, variable=self.graph_props._use_kde)
      self.sep = Checkbutton(self.check_frame, text='Use Separate Plots',
                       onvalue=True, offvalue=False,
                       variable=self.graph_props._separate)
      # Now grid the entries
      self.legend.grid(column=0, row=0, padx=2, pady=2, sticky=N+S+E+W)
      self.gridl.grid(column=1, row=0, padx=2, pady=2, sticky=N+S+E+W)
      self.time.grid(column=2, row=0, padx=2, pady=4, sticky=N+S+E+W)
      self.norm.grid(column=0, row=1, padx=2, pady=4, sticky=N+S+E+W)
      self.kde.grid(column=1, row=1, padx=2, pady=4, sticky=N+S+E+W)
      self.sep.grid(column=2, row=1, padx=2, pady=4, sticky=N+S+E+W)
      self.check_frame.pack(fill=BOTH)
      # Draw the entries
      self.entry_frame = Frame(self)
      self.ti = LabeledEntry(self.entry_frame, 'Plot Title', width=80,
                             textvariable=self.graph_props._title)
      self.xl = LabeledEntry(self.entry_frame, 'X-axis Label', width=80,
                             textvariable=self.graph_props._xlabel)
      self.yl = LabeledEntry(self.entry_frame, 'Y-axis Label', width=80,
                             textvariable=self.graph_props._ylabel)
      self.window = LabeledEntry(self.entry_frame, 'Running Avg. Window',
                       width=20, textvariable=self.graph_props._window)
      self.lw = LabeledEntry(self.entry_frame, 'Line Width', width=20,
                             textvariable=self.graph_props._linewidth)
      self.ps = LabeledEntry(self.entry_frame, 'Point Size', width=20,
                             textvariable=self.graph_props._pointsize)
      self.styles = StyleBoxes(self.entry_frame, self.graph_props._linestyle,
                               self.graph_props._pointstyle)
      self.stride = LabeledEntry(self.entry_frame, 'Stride', width=20,
                                 textvariable=self.graph_props._stride)
      self.bins = LabeledEntry(self.entry_frame, 'Bin Width or -# Bins',
                             width=20, textvariable=self.graph_props._bin_width)
      self.nexcl = LabeledEntry(self.entry_frame, 'Pts Skip From Begin',
                                width=20, textvariable=self.graph_props._nexcl)
      # Now grid the entry buttons
      self.ti.grid(column=0, row=0, columnspan=4, pady=4, sticky=N+S+E+W)
      self.xl.grid(column=0, row=1, columnspan=4, pady=4, sticky=N+S+E+W)
      self.yl.grid(column=0, row=2, columnspan=4, pady=4, sticky=N+S+E+W)
      self.nexcl.grid(column=0, row=3, pady=4, sticky=N+S+E+W)
      self.stride.grid(column=1, row=3, pady=4, sticky=N+S+E+W)
      self.bins.grid(column=2, row=3, pady=4, padx=2, sticky=N+S+E+W)
      self.window.grid(column=3, row=3, pady=4, padx=2, sticky=N+S+E+W)
      self.lw.grid(column=0, row=4, pady=4, padx=2, sticky=N+S+E+W)
      self.ps.grid(column=1, row=4, pady=4, padx=2, sticky=N+S+E+W)
      self.styles.grid(column=2, row=4, pady=4, padx=2, columnspan=2, 
                       sticky=N+S+E+W)
      self.entry_frame.pack(fill=BOTH)

class StyleBoxes(Frame):
   """ Controls the line style and point style giving finite choices """
   
   def __init__(self, master, lsvar, psvar):
      Frame.__init__(self, master)
      # Line style menu
      self.ls = Menubutton(self, text='Line Style', width=20)
      self.lsmenu = Menu(self.ls)
      self.lsmenu.add_radiobutton(label="Solid Line", variable=lsvar, value='-')
      self.lsmenu.add_radiobutton(label="Dashed Line", variable=lsvar,
                                  value='--')
      self.lsmenu.add_radiobutton(label="Dash-Dot Line", variable=lsvar,
                                  value='-.')
      self.lsmenu.add_radiobutton(label="Dotted Line", variable=lsvar,
                                  value=':')
      self.lsmenu.add_radiobutton(label='Different', variable=lsvar,
                                  value='u')
      self.lsmenu.add_radiobutton(label="No line", variable=lsvar, value=' ')
      self.ls.config(menu=self.lsmenu)
      # Point style menu
      self.ps = Menubutton(self, text='Point Style', width=20)
      self.psmenu = Menu(self.ps)
      self.psmenu.add_radiobutton(label="Point", variable=psvar, value='.')
      self.psmenu.add_radiobutton(label="Pixel", variable=psvar, value=',')
      self.psmenu.add_radiobutton(label="Circle", variable=psvar, value='o')
      self.psmenu.add_radiobutton(label="Triangle Down", variable=psvar,
                                  value='v')
      self.psmenu.add_radiobutton(label="Triangle Up", variable=psvar, value='^')
      self.psmenu.add_radiobutton(label="Triangle Left", variable=psvar,
                                  value='<')
      self.psmenu.add_radiobutton(label="Triangle Right", variable=psvar,
                                  value='>')
      self.psmenu.add_radiobutton(label="Square", variable=psvar, value='s')
      self.psmenu.add_radiobutton(label="Star", variable=psvar, value='*')
      self.psmenu.add_radiobutton(label="Pentagon", variable=psvar, value='p')
      self.psmenu.add_radiobutton(label="Hexagon", variable=psvar, value='h')
      self.psmenu.add_radiobutton(label="Hexagon (2)", variable=psvar, value='H')
      self.psmenu.add_radiobutton(label="Plus", variable=psvar, value='+')
      self.psmenu.add_radiobutton(label="Diamond", variable=psvar, value='D')
      self.psmenu.add_radiobutton(label="Different", variable=psvar, value='u')
      self.psmenu.add_radiobutton(label="No Points", variable=psvar, value=' ')
      self.ps.config(menu=self.psmenu)
      # Grid them
      self.ls.grid(column=0, row=0, sticky=N+S+E+W)
      self.ps.grid(column=1, row=0, sticky=N+S+E+W)
