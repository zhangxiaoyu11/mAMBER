"""
Has classes for various windows
"""

from __future__ import division
from Tkinter import *

class TextWindow(Toplevel):
   """ A basic, scrollable text window that can either be editable or not """
   def __init__(self, master, state=DISABLED, noclose=False, safe_close=False,
                width=90, height=30):
      """ Make a scrollable, resizeable text """
      Toplevel.__init__(self, master)
      # Allow disabling of window killing
      if noclose:
         self.protocol('WM_DELETE_WINDOW', self.hideme)
      # See if we want to prompt saving before closing the window
      self.safe_close = safe_close
      # Store our original state so we can restore it after a write()
      self.original_state = state
      # Add a horizontal and vertical scroller
      self.hscroller = Scrollbar(self, orient=HORIZONTAL)
      self.vscroller = Scrollbar(self, orient=VERTICAL)
      # Make the text box
      self.text = Text(self, width=width, state=state, height=height, wrap=NONE,
                       xscrollcommand=self.hscroller.set,
                       yscrollcommand=self.vscroller.set)
      # Pack everything in there, nice and tight. Let the text expand and the
      # scroll bars lengthen, but do not let the scroll bars thicken.
      self.text.grid(column=0, row=0, sticky=N+S+E+W)
      self.hscroller.grid(column=0, row=1, sticky=N+S+E+W)
      self.vscroller.grid(column=1, row=0, sticky=N+S+E+W)
      self.columnconfigure(0, weight=1)
      self.columnconfigure(1, weight=0)
      self.rowconfigure(0, weight=1)
      self.rowconfigure(1, weight=0)
      # Now make the scroll bars actually work
      self.hscroller.configure(command=self.text.xview)
      self.vscroller.configure(command=self.text.yview)

      # Add a 'file' Menu
      filemenu = Menu(self, tearoff=0)
      filemenu.add_command(label='Save', command=self.savetext,
                           accelerator='<Ctrl-s>')
      filemenu.add_command(label='Save As', command=self.saveastext,
                           accelerator='<Ctrl-S>')
      if noclose:
         filemenu.add_separator()
         filemenu.add_command(label='Hide Window', command=self.hideme)
      else:
         filemenu.add_separator()
         filemenu.add_command(label='Close Window', command=self.destroy)
      # Add a menu bar
      menubar = Menu(self, tearoff=0)
      menubar.add_cascade(label='File', underline=0, menu=filemenu)
      self.configure(menu=menubar)

      # Track whether we are hidden or not
      self.hidden = False

      # This is updated as we need to to make sure we know when to ask if we
      # want to save appropriately.  This will only work for a DISABLE'd
      # original state, otherwise we will have to check.
      self.up_to_date = True

      # Bind ctrl-S to the save function
      self.focus_set()
      self.bind("<Control-KeyPress-s>", self.savetext)
      self.bind("<Control-KeyPress-S>", self.saveastext)

   def write(self, s):
      """ 
      Writes 's' to the window, such that it will emulate a file.  We have to
      change the state to ACTIVE in order to add text, but then change it back
      to the original state afterwards
      """
      self.text.configure(state=NORMAL)
      self.text.insert(END, s)
      self.text.configure(state=self.original_state)
      self.up_to_date = False

   def savetext(self, event=None):
      """ Save to the old filename if we already saved once """
      if hasattr(self, 'fname'):
         self.text.configure(state=NORMAL)
         f = open(self.fname, 'w')
         f.write(self.text.get('0.0', END))
         f.close()
         self.text.configure(state=self.original_state)
         self.up_to_date = True

      else:
         return self.saveastext(event)

   def saveastext(self, event=None):
      """ Offer to save this to a new file """
      from tkFileDialog import asksaveasfilename
      from tkMessageBox import showwarning

      fname = asksaveasfilename(parent=self, title='Save Text to File', 
                                defaultextension='.txt')
      # Allow cancelling
      if not fname: return

      self.fname = fname

      try:
         f = open(fname, 'w')
      except IOError:
         showwarning('Bad Permissions', 'Could not open %s for writing!' %
                     fname)
         return
      
      # Activate text object, extract the text, then reset the state
      self.text.configure(state=NORMAL)
      f.write(self.text.get('0.0', END))
      self.text.configure(state=self.original_state)
      f.close()
      self.up_to_date = True

   def destroy(self):
      """ See if we want to save before closing """
      from tkMessageBox import askyesno
      if not self.safe_close: 
         Toplevel.destroy(self)
         return
      if not hasattr(self, 'fname'):
         # No filename?  We have to ask to save
         if askyesno('Save File?', 'Do you want to save the text in this '
                     'window?', parent=self):
            self.saveastext()
         Toplevel.destroy(self)
         return
      # If we are here, we want to ask to save
      if self.original_state is DISABLED and self.up_to_date:
         # Already up-to-date
         Toplevel.destroy(self)
      elif self.original_state is NORMAL or not self.up_to_date:
         # Now we have to do text comparison. Yuck
         file_text = open(self.fname, 'r').read()
         window_text = str(self.text.get('0.0', END))
         if file_text != window_text:
            if askyesno('Save File?', 'Do you want to save the text in this '
                        'window?', parent=self):
               self.savetext()
         Toplevel.destroy(self)
         return

      # See if our text is different
   def hideme(self):
      """ Hides the current window """
      if not self.hidden:
         self.withdraw()
         self.hidden = True

   def revealme(self):
      """ Brings the window back """
      if self.hidden:
         self.deiconify()
         self.hidden = False
