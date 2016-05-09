from Tkinter import *
from mdoutanalyzer.widgets import (DataButton, GraphButton, SaveButton,
                                   StatButton, HistButton, AutoCorrButton,
                                   RunningAvgButton, CumulativeAvgButton)
from mdoutanalyzer.graphproperties import GraphProperties
from mdoutanalyzer.menus import MainMenu

class MdoutAnalyzerApp(Frame):
   """ Main application """
   
   NCOLS  = 4
   NCOLS2 = 3

   def __init__(self, root, mdout):
      """ Set up the main app """
      self.figlist = 0
      self.mdout = mdout
      self.graph_props = GraphProperties()
      Frame.__init__(self, root)
      # Create the main file menu
      self.main_menu = MainMenu(self, self.graph_props, self.mdout)
      root.config(menu=self.main_menu)
      # Create a frame/canvas where we will put all of the Data buttons
      self.data_canv = Frame(self)
      self.data_active = [IntVar() for key in self.mdout.data]
      self.button_list = [DataButton(self.data_canv, self.data_active[i], key)
                          for i, key in enumerate(self.mdout.data)]
      # Label this
      Label(self.data_canv, text='Select Data to Analyze').grid(row=0, column=0,
            columnspan=self.NCOLS, sticky=N+S+E+W)
            
      # Grid the buttons in the canvas
      for i, b in enumerate(self.button_list):
         b.grid(row=1+i//self.NCOLS, column=i%self.NCOLS, sticky=N+S+E+W,
                padx=2, pady=2)
      self.data_canv.grid(column=0, row=0, sticky=N+S+E+W)

      # What do we want to do with the data?
      self.action_canv = Frame(self)
      self.acts = []
      self.acts.append(GraphButton(self.action_canv, self.mdout.data,
                           self.data_active, self.graph_props, 'Graph Them!'))
      self.acts.append(SaveButton(self.action_canv, self.mdout.data,
                           self.data_active, self.graph_props, 'Save to File'))
      self.acts.append(StatButton(self.action_canv, self.mdout.data,
                           self.data_active, self.graph_props, 
                           'Show Statistics'))
      self.acts.append(HistButton(self.action_canv, self.mdout.data,
                           self.data_active, self.graph_props,
                           'Histogram Data'))
      self.acts.append(AutoCorrButton(self.action_canv, self.mdout.data,
                           self.data_active, self.graph_props,
                           'Autocorrelation'))
      self.acts.append(RunningAvgButton(self.action_canv, self.mdout.data,
                           self.data_active, self.graph_props,
                           'Running Average'))
      self.acts.append(CumulativeAvgButton(self.action_canv, self.mdout.data,
                           self.data_active, self.graph_props,
                           'Cumulative Running Average'))

      Label(self.action_canv, 
            text='What do you want to do with your data?').grid(
            row=0, column=0, columnspan=self.NCOLS2, sticky=N+S+E+W)

      # Grid the buttons in the canvas
      for i, b in enumerate(self.acts):
         b.grid(row=1+i//self.NCOLS2, column=i%self.NCOLS2, sticky=N+S+E+W,
                padx=2, pady=2)
      self.action_canv.grid(column=0, row=1, sticky=N+S+E+W)

      # Bind <F1> to the help menu
      self.bind('<F1>', self.Help)
      self.focus_set()

   def Help(self, *args, **kwargs):
      """ Wrapper for the help function """
      Help(self)

   def Credits(self, *args, **kwargs):
      """ Wrapper for the Credits funtion """
      Credits(self)

   def next_fig(self):
      """ Generate the next figure button """
      self.figlist += 1
      return self.figlist

def Help(master):
   """ Displays a help window """

   help_msg = (
      'MdoutAnalyzer.py is a program that parses output files from sander and\n'
      'pmemd and allows you to plot the data from the various energy terms in\n'
      'a number of ways.\n\n'
      
      '  Graph Them! -- Graph the raw data\n'
      '  Save to File -- Save all selected data to separate columns in a text\n'
      '                  or CSV file.\n'
      '  Show Statistics -- Print the min and max values, average, and \n'
      '                     standard deviation of selected data.\n'
      '  Histogram Data -- Plots histograms (optionally smoothed with kernel\n'
      '                    density estimations) of the selectd data\n'
      '  Autocorrelation -- Plots the normalized autocorrelation function of\n'
      '                     the data.\n'
      '  Running Average -- Plots the running average by averaging each point\n'
      '                within a given range on either side of that data point\n'
      '  Cumulative Running Average -- Plots the average of all data points\n'
      '                before and including that data point.\n\n'
      
      'The <Graph Property Control> window provides options to customize the\n'
      'appearance of the plots. This window starts automatically upon startup\n'
      'and can be shown again from the File menu. Some options are described:\n'
      '\n'
      '  Show Legend -- Display the legend on the plot\n'
      '  Show Grid Lines -- Show dotted grid lines at tick marks on plots\n'
      '  Use Time as X-Coord -- Use the TIME(PS) data as the X-coordinate in\n'
      '                         relevant plots.\n'
      '  Normalize Histograms -- Plot a normalized histogram\n'
      '  Smooth Histograms (KDE) -- Use a kernel density estimate to smooth\n'
      '                             histograms. These are ALWAYS normalized.\n'
      '  Use Separate Plots -- Put each selected term in its own subplot\n'
      '  Pts Skip from Begin -- Number of data points to skip from the start\n'
      '                         of the requested data sets.\n'
      '  Stride -- Interval between adjacent plotted points.\n'
      '  Bin Width or -# Bins -- When positive, this is the bin width of the\n'
      '              histograms (or bandwidth of KDEs). When negative, it is\n'
      '              the number of bins. When zero, defaults are used\n'
      '  Running Avg. Window -- # of points on either side of each data point\n'
      '                         to use in computing the running average.\n')

   # Use the message itself to calculate the necessary size for the window
   nlines = len(help_msg.split('\n'))
   maxwidth = max([len(l) for l in help_msg.split('\n')])
   
   # Generate a toplevel window to put the help text into
   window = Toplevel(master)
   window.title('Help')
   text = Text(window, width=maxwidth+1, height=nlines)
   text.insert(END, help_msg)
   text.config(relief=FLAT)
   text.config(state=DISABLED)
   text.pack()
   window.resizable(False, False)
   window.grab_set()

def Credits(master):
   """ Display the credits (About) for this program """
   from mdoutanalyzer.__init__ import __version__, __author__, __date__
   from os import path

   CRW = 450 # credits window width
   CRH = 350 # credits window height

   properties = {'version' : __version__, 'author' : __author__,
                 'date' : __date__}
   header_text = (
      'MdoutAnalyzer.py is a program that parses sander and pmemd output\n'
      'files and allows rapid visualization of output data.\n'
                 )
   tailer_text = (
      'MdoutAnalyzer.py Version %(version)s\n'
      'Written by %(author)s\n'
      ' %(date)s '
                 ) % properties

   credits = Toplevel(master)
   credits.resizable(False, False)
   credits.title('About MdoutAnalyzer.py')
   # Get the file name
   fname = path.join(path.abspath(path.split(__file__)[0]), 'energy2.gif')
   canv = Canvas(credits, width=CRW, height=CRH)
   canv.pack()
   image = PhotoImage(name='SPAM', master=canv, file=fname)
   canv.create_image(CRW//2, CRH//2, anchor=CENTER, image=image)
   canv.create_text(CRW//2, 10, anchor=N, text=header_text,
                    justify=CENTER)
   canv.create_text(CRW//2, CRH-10, anchor=S, text=tailer_text,
                    justify=CENTER)
   canv.update_idletasks()
   credits.grab_set()
   # Save a reference to the image so it doesn't disappear
   master.image = image
   
