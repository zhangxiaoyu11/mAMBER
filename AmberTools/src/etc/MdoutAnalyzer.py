#! PYTHONEXE

from tkFileDialog import askopenfilenames
from Tkinter import Tk, BOTH
import matplotlib
matplotlib.use('TkAgg')
from mdoutanalyzer import __version__, __author__, __date__
from mdoutanalyzer.graphproperties import GraphControlWindow
from mdoutanalyzer.mdout import AmberMdout
from mdoutanalyzer.toplevel_app import MdoutAnalyzerApp
from optparse import OptionParser
import re
import sys

geore = re.compile(r'(\d+)x(\d+)\+(\d+)\+(\d+)')
geoformat = '%dx%d+%d+%d'

verstring = """
   %%prog : An AMBER MD output file parser and graphing utility

                              Version %s
                             %s

   Written by %s
""" % (__version__, __date__, __author__)

parser = OptionParser(usage='%prog [mdout1] [mdout2] ... [mdoutN]',
                      version=verstring)

opt, arg = parser.parse_args()

root = Tk()
root.title('Mdout Analyzer')
if not arg:
   arg = askopenfilenames(title='Select Mdout File(s)', parent=root,
                          filetypes=[('Mdout Files', '*.mdout'),
                                     ('All Files', '*')])

if not arg:
   print ('No mdout files chosen!')
   sys.exit(1)

for f in arg:
   try:
      mdout += AmberMdout(f)
   except NameError:
      mdout = AmberMdout(f)

app = MdoutAnalyzerApp(root, mdout)
app.pack(fill=BOTH)
# Update idle tasks here to make sure the whole app is filled in before making
# our window non-resizable.  In some instances, this can chop off the second
# frame
app.update_idletasks()
# Now make our window non-resizable
root.resizable(False, False)
root.update_idletasks()
# Load up the graph options window and move it to the right of the root window
rootgeo = [int(i) for i in geore.match(root.geometry()).groups()]
graphprops = GraphControlWindow(root, app.graph_props)
root.update_idletasks()
ggeo = [int(i) for i in geore.match(graphprops.geometry()).groups()]
graphprops.geometry(geoformat % (ggeo[0], ggeo[1], rootgeo[2] + rootgeo[0] + 5,
                                 rootgeo[3])
                   )

# For some reason it seems like mainloop is not exited properly if we simply
# call root.destroy(). Therefore, we'll replace WM_DELTE_WINDOW with root.quit()
# to bust out of the mainloop.
root.protocol('WM_DELETE_WINDOW', root.quit)

# Enter our mainloop
root.mainloop()
