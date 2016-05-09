from Tkinter import *
from tkFileDialog import askopenfilenames
from mdout import AmberMdout
from mdoutanalyzer.graphproperties import GraphControlWindow

class FileMenu(Menu):
   """ The main file menu """
   def __init__(self, master, mdout, graph_props):
      Menu.__init__(self, master, tearoff=0)
      self.mdout = mdout
      self.graph_props = graph_props
      self.add_command(label='Add Mdout File(s)', command=self._add_mdout)
      self.add_separator()
      self.add_command(label='Edit Graph Properties', command=self._edit_props)
      self.add_separator()
      self.add_command(label='Quit', command=master.master.quit)

   def _add_mdout(self):
      """ Get open filenames """
      fnames = askopenfilenames(title='Select Mdout File(s)', parent=self,
                                filetypes=[('Mdout Files', '*.mdout'),
                                           ('All Files', '*')])
      for f in fnames:
         self.mdout += AmberMdout(f)

   def _edit_props(self):
      """ Open up a window controlling the graph properties """
      GraphControlWindow(self, self.graph_props)

class HelpMenu(Menu):
   """ About and Help """
   def __init__(self, master):
      Menu.__init__(self, master, tearoff=0)
      # Holders
      self.add_command(label='About', command=master.master.Credits)
      self.add_separator()
      self.add_command(label='Help', command=master.master.Help)

class MainMenu(Menu):
   """ The main file menu for the main App """

   def __init__(self, master, graph_props, mdout):
      """ Sets up the main menu """
      Menu.__init__(self, master, tearoff=0, relief=FLAT,
                    activeborderwidth=0)
      self.graph_props = graph_props
      self.mdout = mdout
      self.file_menu = FileMenu(self, self.mdout, self.graph_props)
      self.help = HelpMenu(self)
      # Add the menus
      self.add_cascade(label='File', underline=0, menu=self.file_menu)
      self.add_cascade(label='Help', underline=0, menu=self.help)

if __name__ == '__main__':
   from mdoutanalyzer.graphproperties import GraphProperties
   # Testing part

   root = Tk()
   root.config(menu=MainMenu(root, GraphProperties(), object()))
   root.mainloop()
