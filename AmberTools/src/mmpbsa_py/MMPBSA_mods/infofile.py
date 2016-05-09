"""
This module contains classes that back up the state of the MM/PBSA calculation
so that future post-processing of the output files can be done without
re-supplying all of the information again.

######################### GPL LICENSE INFO ############################

  Copyright (C) 2009 - 2012  Jason Swails, Bill Miller III, and Dwight McGee

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
   
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330,
  Boston, MA 02111-1307, USA.

"""

import re
import warnings

class InfoFile(object):
   """
   This is the class responsible for backing up all of the calculation
   information and details so the output files can be re-written or the data
   can be readied for post-processing with external scripts.

   This object is passed a copy of the main app MMPBSA_App class
   """
   
   # Make a list about which INPUT variables are editable
   EDITABLE_INFO_VARS = ['debug_printlevel', 'verbose', 
                         'csv_format', 'dec_verbose']

   def __init__(self, app):
      """ Instantiate me """
      self.app = app

   def write_info(self, name=None):
      """ Writes the info file to the info name """
      if name is None:
         name = self.app.pre + 'info'
      outfile = open(name, 'w')
      # The data we have to write: INPUT, FILES, and the following attributes:
      # numframes, numframes_nmode, mpi_size (also recognize size), 
      # input_file_text, mut_str

      # Start with INPUT (and the editable vars). Allow this to recognize INFO
      # files from the last version of MMPBSA.py
      outfile.write('# You can alter the variables below to change what info '
                    'is printed out\n')
      for var in InfoFile.EDITABLE_INFO_VARS:
         outfile.write("INPUT['%s'] = %s\n" % (var, 
               self.write_var(self.app.INPUT[var])))
      outfile.write('#\n# MODIFY NOTHING BELOW HERE, OR GET WHAT YOU DESERVE\n')
      for var in self.app.INPUT.keys():
         outfile.write("INPUT['%s'] = %s\n" % (var,
               self.write_var(self.app.INPUT[var])))

      # Now print out the FILES
      for var in dir(self.app.FILES):
         # Skip over __method__ functions and output files
         if var.startswith('__') or var in ('rewrite_output', 'output_file', 
                  'decompout', 'energyout', 'dec_energies', 'overwrite'):
            continue
         outfile.write("FILES.%s = %s\n" % (var, 
                                  self.write_var(getattr(self.app.FILES, var))))

      # Now print out the attributes we need to save
      outfile.write('size = %d\n' % self.app.mpi_size)
      outfile.write('numframes = %d\n' % self.app.numframes)
      outfile.write('numframes_nmode = %d\n' % self.app.numframes_nmode)
      outfile.write("mut_str = '%s'\n" % self.app.mut_str)
      outfile.write(self.app.input_file_text)

   def read_info(self, name=None):
      """ Reads a _MMPBSA_info file and populates the app """
      from MMPBSA_mods.commandlineparser import OptionList
      # Give self.app a FILES attribute if it does not have one yet before
      # trying to modify it
      if not hasattr(self.app, 'FILES'):
         self.app.FILES = OptionList()
      if name is None:
         name = self.app.pre + 'info'
      #TODO Add checking for differences in FILES
      inputre = re.compile(r'''INPUT\['(\S+)'\] = ((?:'.*')|(?:".*")|\S+)''')
      filesre = re.compile(r'''FILES\.(\S+) = (.*)''')
      otherre = re.compile(r'(\S+) = (\S+)')
      infile = open(name, 'r')
      input_text = ''
      for line in infile:
         if line.startswith('#'):
            continue
         # Now load each info file line
         if line.startswith('|'):
            input_text += line
            continue
         rematch = inputre.match(line)
         if rematch:
            var, val = rematch.groups()
            val = _determine_type(val)
            self.app.INPUT[var] = val
            continue
         rematch = filesre.match(line)
         if rematch:
            var, val = rematch.groups()
            val = _determine_type(val)
            setattr(self.app.FILES, var, val)
         rematch = otherre.match(line)
         if rematch:
            var, val = rematch.groups()
            val = _determine_type(val)
            if var == 'size':
               self.app.mpi_size = val
            else:
               setattr(self.app, var, val)
            continue
      # Determine stability here:
      self.app.stability = (self.app.FILES.receptor_prmtop is None and
                            self.app.FILES.ligand_prmtop is None)
      if self.app.FILES.receptor_mdcrd or self.app.FILES.ligand_mdcrd:
         self.app.traj_protocol = 'MTP'
      else:
         self.app.traj_protocol = 'STP'
      # Load the input file text
      self.app.input_file_text = input_text

   def write_var(self, var):
      """ 
      Wrapper to return a string in which str vars are enclosed in quotes and
      numeric types (int and float) are not
      """
      if isinstance(var, str):
         return "'%s'" % var
      else:
         return "%s" % var


def _determine_type(thing):
   """ Determines what type this thing is """
   # If it is a string with quotes, strip off the quotes
   if thing.startswith('"') and thing.endswith('"'):
      return thing[1:len(thing)-1]
   if thing.startswith("'") and thing.endswith("'"):
      return thing[1:len(thing)-1]

   # Check for bool or None
   if thing == 'True':
      return True
   elif thing == 'False':
      return False
   elif thing == 'None':
      return None

   # Check for int, then check for float
   try:
      return int(thing)
   except ValueError:
      pass

   try:
      return float(thing)
   except ValueError:
      pass

   # Check for list
   if thing.startswith('[') and thing.endswith(']'):
      return eval(thing)

   # No idea what else it could be! Return string, but warn
   warnings.warn('Encountered unknown type in info file.\n' +
                 str(thing))
   return thing
