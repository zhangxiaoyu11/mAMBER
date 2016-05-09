"""
List of topology file objects that can be edited in ParmEd. They can be indexed
via either the name of the topology file or by the order in which they were
loaded.
"""

from chemistry.amber.readparm import AmberParm
from ParmedTools.exceptions import DuplicateParm

class ParmList(object):
   """
   List of AmberParm objects index-able via either prmtop name or prmtop index
   (based on order added)
   """
   def __init__(self):
      """ Must be instantiated by itself """
      self._parm_names = []
      self._parm_instances = []
      self.parm = None # The active parm instance
      self.current_index = 0

   def set_new_active(self, idx):
      """ Sets a new active parm """
      self.current_index = self.index(idx)
      self.parm = self[self.current_index]

   def add_parm(self, parm):
      """ Add a parm to the list """
      # Convert a string to an AmberParm or add an AmberParm directly
      if isinstance(parm, str):
         parm = AmberParm(parm)
      elif not isinstance(parm, AmberParm):
         raise TypeError('ParmList can only add AmberParm instances')
      # Make sure this parm is not part of the list already
      if str(parm) in self._parm_names:
         raise DuplicateParm('%s already in ParmList' % parm)
      # Otherwise, add in the new parm's name
      self._parm_names.append(str(parm))
      self._parm_instances.append(parm)
      # A newly added topology file is the currently active parm
      self.current_index = len(self._parm_instances) - 1
      self.parm = parm

   def index(self, idx):
      """ See what the index of the requested parm is (can be int or str) """
      if isinstance(idx, int):
         if idx < 0 or idx >= len(self._parm_instances):
            raise IndexError('index out of range for ParmList')
         return idx
      elif isinstance(idx, str):
         try:
            return self._parm_names.index(idx)
         except ValueError:
            raise IndexError('%s prmtop not in ParmList' % idx)

      raise TypeError('ParmList is only indexed via parm names and integers')

   def __getitem__(self, idx):
      """ Allow an integer index or string index to identify a parm """
      return self._parm_instances[self.index(idx)]

   def __contains__(self, name):
      """ See if a given parm name is in the list """
      if isinstance(name, str):
         return name in self._parm_names
      elif isinstance(name, int):
         return name >= 0 and name < len(self)
      else:
         return False
   
   def __len__(self):
      return len(self._parm_instances)

   def empty(self):
      return len(self._parm_instances) == 0
