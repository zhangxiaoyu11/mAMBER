"""
Module for evaluating Amber Mask strings and translating them into lists in which
a selected atom is 1 and one that's not is 0.
"""
from chemistry.exceptions import MaskError
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class AmberMask(object):
   """ 
   What is hopefully a fully-fledged Amber mask parser implemented in Python.

   It currently lacks the capability to evaluate masks with distance criteria
   """

   #======================================================

   def __init__(self, parm, mask):
      self.parm = parm
      self.mask = mask
      
   #======================================================

   def __str__(self):
      return self.mask

   #======================================================

   def Selected(self, invert=False):
      """ Generator that returns the indexes of selected atoms """
      for i, v in enumerate(self.Selection(invert=invert)):
         if v:
            yield i

   #======================================================

   def Selection(self, prnlev=0, invert=False):
      """ Parses the mask and analyzes the result to return an atom
          selection array
      """
      from sys import stderr, stdout
      if prnlev > 2: stderr.write('In AmberMask.Selection(), debug active!\n')
      if prnlev > 5: stdout.write('original mask: ==%s==\n' % self.mask)

      # 0) See if we got the default "everything" mask(*) and return accordingly
      if self.mask.strip() == '*':
         return [1 for i in range(self.parm.ptr('natom'))]

      # 1) preprocess input expression
      infix = AmberMask._tokenize(self, prnlev)
      if prnlev > 5: stdout.write('tokenized mask: ==%s==\n' % infix)

      # 2) construct postfix (RPN) notation
      postfix = AmberMask._torpn(self, infix, prnlev)
      if prnlev > 5: stdout.write('postfix mask: ==%s==\n' % postfix)

      # 3) evaluate the postfix notation
      if invert:
         return [abs(i-1) for i in AmberMask._evaluate(self, postfix, prnlev)]
      return AmberMask._evaluate(self, postfix, prnlev)

   #======================================================

   def _tokenize(self, prnlev):
      """ Tokenizes the mask string into individual selections:
          1. remove spaces
          2. isolate 'operands' into brackets [...]
          3. split expressions of the type :1-10@CA,CB into 2 parts;
             the 2 parts are joned with & operator and (for the sake
             of preserving precedence of other operators) enclosed by
             (...); i.e. :1-10@CA,CB is split into (:1-10 & @CA,CB)
          4. do basic error checking
      """
      buffer = '' # keeping track of a single operand
      infix = ''  # value that is returned at the end
      # flag == 0: means new operand or operand was completed & ended with ]
      # flag == 1: means operand with ":" read
      # flag == 2: means operand with "@" read
      # flag == 3: means '<' or '>' read, waiting for numbers
      flag = 0
      i = 0
      while i < len(self.mask):
         p = self.mask[i]
         # skip whitespace
         if p.isspace(): 
            i += 1
            continue
         # If p is an operator, is the last character, or is a ()...
         elif self._isOperator(p) or i == len(self.mask) - 1 or p in ['(',')']:
            # Deal with the last character being a wildcard that we have to
            # convert
            if p == '=' and i == len(self.mask) - 1: # wildcard
               if flag > 0: p = '*'
               else: raise MaskError('AmberMask: \'=\' not in name list syntax')
            # If this is the end of an operand, terminate the buffer, flush
            # it to infix, and reset flag to 0 and empty the buffer
            if flag > 0:
               if i == len(self.mask) - 1 and p != ')': buffer += p
               buffer += '])'
               flag = 0
               infix += buffer
               buffer = ''
            if i != len(self.mask) - 1 or p == ')': infix += p
            # else if p is >,<
            if p in ['<','>']:
               buffer = '([%s' % p
               i += 1
               p = self.mask[i]
               buffer += p
               flag = 3
               try:
                  self.parm.coords[0]
               except AttributeError:
                  raise MaskError('AmberMask: <,> operators require ' +
                        'a loaded restart file!')
               if not p in [':','@']:
                  raise MaskError('AmberMask: Bad syntax')
         elif self._isOperand(p):
            if flag == 0:
               buffer = '(['
               flag = 1
               if p != '*':
                  raise MaskError('AmberMask: Bad syntax')
            if p == '=': # wildcard
               if flag > 0: p = '*'
               else: raise MaskError('AmberMask: \'=\' not in name list syntax')
            buffer += p
         elif p == ':':
            if flag == 0:
               buffer = '([:'
               flag = 1
            else:
               buffer += '])|([:'
               flag = 1
         elif p == '@':
            if flag == 0:
               buffer = '([@'
               flag = 2
            elif flag == 1:
               buffer += ']&[@'
               flag = 2
            elif flag == 2:
               buffer += '])|([@'
               flag = 2
         else:
            raise MaskError('AmberMask: Unknown symbol (%s) expression when parsing Mask' % p)
         i += 1
      # end while i < len(self.mask):
      # Check that each operand has at least 4 characters: [:1] and [@C], etc.
      i = 0
      n = 1 # number of characters in current operand
      flag = 0
      while i < len(infix):
         p = infix[i]
         if p == '[':
            n += 1
            flag = 1
         elif p == ']':
            if n < 4 and infix[i-1] != '*':
               raise MaskError('AmberMask: empty token in infix')
            n = 1
         else:
            if flag == 1:
               n += 1
         i += 1
      
      return infix + '_' # terminating _ for next step

   #======================================================
   
   def _isOperator(self, char):
      """ Determines if a character is an operator """
      return len(char) == 1 and char in '!&|<>'

   #======================================================

   def _isOperand(self, char):
      """ Determines if a character is an operand """
      return len(char) == 1 and (char in "*/%-?,'.=+" or char.isalnum())
   
   #======================================================

   def _torpn(self, infix, prnlev):
      """ Converts the infix to an RPN array """
      postfix = ''
      stack = ['_']  # use a list as a stack. Then pop() works as expected
      flag = 0
      i = 0

      while i < len(infix):
         p = infix[i]
         if p == '[':
            postfix += p
            flag = 1
         elif p == ']':
            postfix += p
            flag = 0
         elif flag: postfix += p
         elif p == '(': stack.append(p)
         elif p == ')':
            pp = stack.pop()
            while pp != '(':
               if pp == '_':
                  raise MaskError('AmberMask: Unbalanced parentheses in Mask.')
               postfix += pp
               pp = stack.pop()
         # At this point both ()s are discarded
         elif p == '_':
            pp = stack.pop()
            while pp != '_':
               if pp == '(':
                  raise MaskError('AmberMask: Unbalanced parentheses in Mask.')
               postfix += pp
               pp = stack.pop()
         elif self._isOperator(p):
            P1 = self._priority(p)
            P2 = self._priority(stack[len(stack)-1])
            if P1 > P2:
               stack.append(p)
            else:
               while P1 <= P2:
                  pp = stack.pop()
                  postfix += pp
                  P1 = self._priority(p)
                  P2 = self._priority(stack[len(stack)-1])
               stack.append(p)
         else:
            raise MaskError('AmberMask: Unknown symbol %s' % p)
         i += 1
      # end while i < len(infix):
      return postfix

   #======================================================

   def _evaluate(self, postfix, prnlev):
      """ Evaluates a postfix in RPN format and returns a selection array """
      from sys import stderr
      buffer = ''
      stack = []

      pos = 0 # position in postfix
      while pos < len(postfix):
         p = postfix[pos]
         if p == '[': buffer = ''
         elif p == ']': # end of the token
            ptoken = buffer
            pmask = self._selectElemMask(ptoken)
            stack.append(pmask)
         elif self._isOperand(p) or p in [':','@']:
            buffer += p
         elif p in ['&','|']:
            pmask1 = None
            pmask2 = None
            try:
               pmask1 = stack.pop()
               pmask2 = stack.pop()
               pmask = self._binop(p, pmask1, pmask2)
            except IndexError:
               raise MaskError('AmberMask: Illegal binary operation')
            stack.append(pmask)
         elif p in ['<','>']:
            if pos < len(postfix)-1 and postfix[pos+1] in [':','@']: buffer += p
            else:
               try:
                  pmask1 = stack.pop() # this should be the distance criteria
                  pmask2 = stack.pop()
                  pmask = self._selectDistd(pmask1, pmask2)
               except IndexError:
                  return [0 for i in range(self.parm.ptr('natom'))]
               stack.append(pmask)
         elif p == '!':
            try: pmask1 = stack.pop()
            except IndexError: raise MaskError('AmberMask: Illegal ! operation')
            pmask = self._neg(pmask1)
            stack.append(pmask)
         else:
            raise MaskError('AmberMask: Unknown symbol evaluating RPN: %s' % p)
         pos += 1
      # end while i < len(postfix)

      pmask = stack.pop()

      if stack: raise MaskError('AmberMask: There may be missing operands in the mask!')

      if prnlev > 7: stderr.write('%d atoms selected by %s' % (sum(pmask), self.mask))

      return pmask

   #======================================================
   
   def _neg(self, pmask1):
      """ Negates a given mask """
      return pmask1.Not()

   #======================================================

   def _selectDistd(self, pmask1, pmask2):
      """ Selects atoms based on a distance criteria """
      raise MaskError("Distance mask is not implemented yet!")
      pass

   #======================================================

   def _selectElemMask(self, ptoken):
      """ Selects an element mask """
      # some constants
      ALL = 0
      NUMLIST = 1
      NAMELIST = 2
      TYPELIST = 3
      ELEMLIST = 4
      # define the mask object and empty buffer
      pmask = _mask(self.parm.ptr('natom'))
      buffer = ''
      buffer_p = 0
      # This is a residue NUMber LIST
      if ptoken.startswith(':'):
         reslist = NUMLIST
         pos = 1
         while pos < len(ptoken):
            p = ptoken[pos]
            buffer += p
            buffer_p += 1
            if p == '*':
               if buffer_p == 1 and (pos == len(ptoken)-1 or ptoken[pos+1] == ','):
                  reslist = ALL
               elif reslist == NUMLIST: reslist = NAMELIST
            elif p.isalpha() or p == '?':
               reslist = NAMELIST
            if pos == len(ptoken) - 1: buffer_p = 0
            if len(buffer) != 0 and buffer_p == 0:
               if reslist == ALL: pmask.select_all()
               elif reslist == NUMLIST: self._residue_numlist(buffer, pmask)
               elif reslist == NAMELIST: self._residue_namelist(buffer, pmask)
               reslist = NUMLIST
            pos += 1
      elif ptoken.startswith('@'):
         atomlist = NUMLIST
         pos = 1
         while pos < len(ptoken):
            p = ptoken[pos]
            buffer += p
            buffer_p += 1
            if p == '*':
               if buffer_p == 0 and (pos == len(ptoken)-1 or ptoken[pos+1] == ','):
                  atomlist = ALL
               elif atomlist == NUMLIST: atomlist = NAMELIST
            elif p.isalpha() or p == '?':
               if atomlist == NUMLIST: atomlist = NAMELIST
            elif p == '%': atomlist = TYPELIST
            elif p == '/': atomlist = ELEMLIST
            if pos == len(ptoken) - 1: buffer_p = 0

            if len(buffer) != 0 and buffer_p == 0:
               if atomlist == ALL: pmask.select_all()
               elif atomlist == NUMLIST: self._atom_numlist(buffer, pmask)
               elif atomlist == NAMELIST: self._atom_namelist(buffer, pmask)
               elif atomlist == TYPELIST: self._atom_typelist(buffer[1:], pmask)
               elif atomlist == ELEMLIST: self._atom_elemlist(buffer[1:], pmask)
            pos += 1
      elif ptoken.strip() == '*':
         pmask.select_all()
      elif ptoken[0] in ['<','>']:
         return _mask(self.parm.ptr('natom')) # empty mask, this is ignored anyway
      else: raise MaskError('AmberMask: Mask is missing : and @')
      # end if ':' in ptoken:

      return pmask

   #======================================================

   def _atom_numlist(self, instring, mask):
      """ Fills a _mask based on atom numbers """
      buffer = ''
      pos = 0
      at1 = at2 = dash = 0
      while pos < len(instring):
         p = instring[pos]
         if p.isdigit():
            buffer += p
         if p == ',' or pos == len(instring) - 1:
            if dash == 0: 
               at1 = int(buffer)
               self._atnum_select(at1, at1, mask)
            else:
               at2 = int(buffer)
               self._atnum_select(at1, at2, mask)
               dash = 0
            buffer = ''
         elif p == '-':
            at1 = int(buffer)
            dash = 1
            buffer = ''
         if not (p.isdigit() or p in [',','-']):
            raise MaskError('AmberMask: Unknown symbol in atom number parsing')
         pos += 1

   #======================================================

   def _atom_namelist(self, instring, mask, key='ATOM_NAME'):
      """ Fills a _mask based on atom names/types """
      buffer = ''
      pos = 0
      while pos < len(instring):
         p = instring[pos]
         if p.isalnum() or p in ['*','?','+',"'",'-']:
            buffer += p
         if p == ',' or pos == len(instring) - 1:
            if '-' in buffer and buffer[0].isdigit():
               self._atom_numlist(buffer, mask)
            else:
               self._atname_select(buffer, mask, key)
            buffer = ''
         if not (p.isalnum() or p in ",?*'+-"):
            raise MaskError('AmberMask: Unrecognized symbol in atom name parsing %s' % p)
         pos += 1

   #======================================================

   def _atom_typelist(self, buffer, mask):
      """ Fills a _mask based on atom types """
      self._atom_namelist(buffer, mask, key='AMBER_ATOM_TYPE')

   #======================================================

   def _atom_elemlist(self, buffer, mask):
      """ Fills a _mask based on atom elements. For now it will just be Atom
          names, since elements are not stored in the prmtop anywhere.
      """
      self._atom_namelist(self, buffer, mask, key='ATOM_NAME')

   #======================================================

   def _residue_numlist(self, instring, mask):
      """ Fills a _mask based on residue numbers """
      buffer = ''
      pos = 0
      at1 = at2 = dash = 0
      while pos < len(instring):
         p = instring[pos]
         if p.isdigit():
            buffer += p
         if p == ',' or pos == len(instring) - 1:
            if dash == 0: 
               at1 = int(buffer)
               self._resnum_select(at1, at1, mask)
            else:
               at2 = int(buffer)
               self._resnum_select(at1, at2, mask)
               dash = 0
            buffer = ''
         elif p == '-':
            at1 = int(buffer)
            dash = 1
            buffer = ''
         if not (p.isdigit() or p in [',','-']):
            raise MaskError('AmberMask: Unknown symbol in residue number parsing')
         pos += 1

   #======================================================

   def _residue_namelist(self, instring, mask):
      """ Fills a _mask based on residue names """
      buffer = ''
      pos = 0
      while pos < len(instring):
         p = instring[pos]
         if p.isalnum() or p in ['*','?','+',"'",'-']:
            buffer += p
         if p == ',' or pos == len(instring) - 1:
            if '-' in buffer and buffer[0].isdigit():
               self._residue_numlist(buffer, mask)
            else:
               self._resname_select(buffer, mask)
            buffer = ''
         if not (p.isalnum() or p in ",?*'+-"):
            raise MaskError('AmberMask: Unrecognized symbol in residue name parsing')
         pos += 1
      

   #======================================================
   
   def _atnum_select(self, at1, at2, mask):
      """ Fills a _mask array between atom numbers at1 and at2 """
      for i in range(at1-1, at2): mask[i] = 1

   #======================================================
   
   def _resnum_select(self, res1, res2, mask):
      """ Fills a _mask array between residues res1 and res2 """
      for i in range(self.parm.ptr('natom')):
         res = self.parm.residue_container[i]
         if res >= res1 and res <= res2: mask[i] = 1

   #======================================================
   
   def _atname_select(self, atname, mask, key='ATOM_NAME'):
      """ Fills a _mask array with all atom names of a given name """
      for i in range(self.parm.ptr('natom')):
         if _nameMatch(atname, self.parm.parm_data[key][i]):
            mask[i] = 1
         elif atname.isdigit():
            mask[i] = int(int(atname) == i+1)

   #======================================================
   
   def _resname_select(self, resname, mask):
      """ Fills a _mask array with all residue names of a given name """
      for i in range(self.parm.ptr('natom')):
         if _nameMatch(resname, self.parm.parm_data['RESIDUE_LABEL'][
                       self.parm.residue_container[i]-1]):
            mask[i] = 1
         elif resname.isdigit():
            mask[i] = int(int(resname) == self.parm.residue_container[i])
            
   #======================================================
   
   def _binop(self, op, pmask1, pmask2):
      """ Does a binary operation on a pair of masks """
      if op == '&': return pmask1.And(pmask2)
      if op == '|': return pmask1.Or(pmask2)
      raise MaskError('AmberMask: Unknown operator (%s)' % op)

   #======================================================

   def _priority(self, op):
      if op in ['>','<']: return 6
      if op in ['!']: return 5
      if op in ['&']: return 4
      if op in ['|']: return 3
      if op in ['(']: return 2
      if op in ['_']: return 1

      raise MaskError('AmberMask: Unknown operator (%s) in Mask ==%s==' % (op, self.mask))

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def _nameMatch(atnam1, atnam2):
   """ Determines if atnam1 matches atnam2, where atnam1 can have * as a wildcard
       and spaces are ignored. atnam2 should come from the prmtop. We'll use regex
       to do this. 
       
       We will replace * with a regex that will match any alphanumeric character
       0 or more times: * --> \\w*

       We will replace ? with a regex that will match exactly 1 alphanumeric
       character: ? --> \\w

       Then, we will substitute all instances of atnam2 in atnam1 with ''. If it's
       a complete match, then our result will be a blank string (and will evaluate
       to False for boolean conditions). If it's not blank, then it's not a full
       match and should return False
   """
   import re
   atnam1 = str(atnam1).replace(' ','')
   atnam2 = str(atnam2).replace(' ','')
   # Replace amber mask wildcards with appropriate regex wildcards and protect the +
   atnam1 = atnam1.replace('*',r'\S*')
   atnam1 = atnam1.replace('?',r'\S')
   atnam1 = atnam1.replace('+',r'\+')
   # Now replace just the first instance of atnam2 in atnam2 with '', and return *not* that
   # DEBUG:
#  print 'Comparing ==%s== with ==%s==' % (atnam1, atnam2)
#  print not bool(re.sub(atnam1, '', atnam2, 1))
   # END DEBUG
   return not bool(re.sub(atnam1, '', atnam2, 1))

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class _mask(list):
   """ Mask array; only used by AmberMask """

   def __init__(self, natom):
      self.natom = natom
      list.__init__(self, [0 for i in range(natom)])

   def append(self, *args, **kwargs):
      raise MaskError('_mask is a fixed-length array!')

   def extend(self, *args, **kwargs):
      raise MaskError('_mask is a fixed-length array!')

   def pop(self, *args, **kwargs):
      return self[len(self)-1]

   def remove(self, *args, **kwargs):
      raise MaskError('_mask is a fixed-length array!')

   def And(self, other):
      if self.natom != other.natom: 
         raise MaskError("_mask: and() requires another mask of equal size!")
      new_mask = _mask(self.natom)
      for i in range(len(self)):
         new_mask[i] = int(self[i] and other[i])
      return new_mask

   def Or(self, other):
      if self.natom != other.natom:
         raise MaskError('_mask: or() requires another mask of equal size!')
      new_mask = _mask(self.natom)
      for i in range(len(self)):
         new_mask[i] = int(self[i] or other[i])
      return new_mask
   
   def Not(self):
      new_mask = _mask(self.natom)
      for i in range(self.natom):
         new_mask[i] = int(not self[i])
      return new_mask

   def select_all(self):
      for i in range(self.natom):
         self[i] = 1

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
