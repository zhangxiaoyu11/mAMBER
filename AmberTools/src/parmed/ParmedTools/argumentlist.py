"""
This stores a list of arguments, tokenizing a string into a list of arguments.
"""
import warnings
from ParmedTools.exceptions import NoArgument, ParmError

class Argument(object):
   """ Is an argument """
   def __init__(self, string):
      """ Input string for the argument """
      self.string = string.strip()
      self.marked = False
      self.keep_quotes = False
      if len(string) == 0 or (len(string) == 1 and string in ("'", '"')):
         warnings.warn("Unclosed quotation detected! Ignoring lone quote")

   def lower(self):
      return self.string.lower()

   def isfloat(self):
      """ Determines if it can be a double """
      try:
         tmp = float(self.string)
         del tmp
         return True
      except ValueError:
         return False

   def isint(self):
      """ Determines if it can be an integer """
      try:
         tmp = int(self.string)
         del tmp
         return True
      except ValueError:
         return False

   def ismask(self):
      """ Determines if we can be an Amber mask or not """
      # The easiest way to do this is to just search for ":" or "@" in the
      # input string (or just the string *, which matches everything).
      return self.string == '*' or ':' in self.string or '@' in self.string

   # Define some casting behavior
   def __float__(self):
      self.marked = True
      return float(self.string)

   def __int__(self):
      self.marked = True
      return int(self.string)

   def __str__(self):
      self.marked = True
      if self.keep_quotes:
         return self.string

      # Strip leading and trailing quotes (only if they both exist)
      if self.string[0] == "'" and self.string[len(self.string)-1] == "'":
         return self.string[1:len(self.string)-1]
      elif self.string[0] == '"' and self.string[len(self.string)-1] == '"':
         return self.string[1:len(self.string)-1]

      return self.string

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class ArgumentList(object):
   """ 
   List of arguments.  
   The arguments should be parsed in the following order to ensure it's done
   correctly:
      get_key_*
      get_next_int
      get_next_float
      get_next_mask
      get_next_string
   """
   def __init__(self, inp):
      """ Defines the string """
      if isinstance(inp, str) or isinstance(inp, int) or isinstance(inp, float):
         # pad with spaces for regex and quotation marks
         self.original = ' %s ' % inp
         self.tokenize(self.original)
      elif isinstance(inp, tuple) or isinstance(inp, list):
         inp = [l for l in inp if l is not None and l != '']
         self.original = ' '.join([str(l) for l in inp])
         self.tokenlist = [Argument(str(l)) for l in inp]
      else:
         raise TypeError("Unrecognized type for ArgumentList: %s" % type(inp))

   def tokenize(self, instring):
      """ 
      Tokenizes the line, everything between quotes is a single token,
      including whitespace
      """
      import re
      tokenre = re.compile(r'''((?:\s+".*"\s+)|(?:\s+'.*'\s+)|(?:\S+))''')
      tokenlist = tokenre.findall(instring)
      # Make sure every non-whitespace character is consumed
      if len(tokenre.sub('', instring).strip()) > 0:
         warnings.warn("Not all of argument string were handled: [%s]" % 
               tokenre.sub('', instring).strip())
      # Strip out whitespace in all arguments
      self.tokenlist = [Argument(token.strip()) for token in tokenlist]

   def __str__(self):
      """ Return the original input string """
      return self.original

   def get_next_int(self, optional=False, default=None):
      """ Returns the next unmarked integer in the token list """
      for token in self.tokenlist:
         if token.isint() and not token.marked:
            return int(token)
      if optional: return default
      raise NoArgument("Missing integer argument!")

   def get_next_float(self, optional=False, default=None):
      """ Returns the next unmarked float in the token list """
      for token in self.tokenlist:
         if token.isfloat() and not token.marked:
            return float(token)
      if optional: return default
      raise NoArgument("Missing floating point argument!")

   def get_next_mask(self, optional=False, default=None):
      """ Returns the next string as long as it matches a mask """
      for token in self.tokenlist:
         if token.marked: continue
         if token.ismask():
            return str(token)
      if optional: return default
      raise NoArgument("Missing Amber mask!")

   def get_next_string(self, optional=False, keep_quotes=False, default=None):
      """ Returns the next unmarked float in the token list """
      for token in self.tokenlist:
         if not token.marked:
            token.keep_quotes = keep_quotes
            return str(token)
      if optional: return default
      raise NoArgument("Missing string!")

   def unmarked(self):
      """ Returns all unmarked arguments as a list of strings """
      unmarked_tokens = []
      for token in self.tokenlist:
         if not token.marked:
            unmarked_tokens.append(token.string)
      return unmarked_tokens

   def _get_key_arg(self, key, argtype):
      """ Gets a key with a given argument type """
      for i, token in enumerate(self.tokenlist):
         if token.lower() == key.lower():
            if self.tokenlist[i+1].marked:
               raise ParmError("Expected token already marked! Should not " +
                               "happen. This means a buggy action...")
            token.marked = True
            return argtype(self.tokenlist[i+1])
      raise NoArgument('Catch me!')

   def get_key_int(self, key, default):
      """ Get an integer that follows a keyword """
      try:
         return self._get_key_arg(key, int)
      except NoArgument:
         return default

   def get_key_float(self, key, default):
      """ Get a float that follows a keyword """
      try:
         return self._get_key_arg(key, float)
      except NoArgument:
         return default

   def get_key_mask(self, key, default):
      """ Get a mask that follows a keyword """
      for i, token in enumerate(self.tokenlist):
         if token.lower() == key.lower():
            if not self.tokenlist[i+1].ismask():
               raise TypeError("Expected mask to follow %s. Got %s instead" % (
                               key, self.tokenlist[i+1].string))
      try:
         return self._get_key_arg(key, str)
      except NoArgument:
         return default

   def get_key_string(self, key, default):
      """ Get a string that follows a keyword """
      try:
         return self._get_key_arg(key, str)
      except NoArgument:
         return default
