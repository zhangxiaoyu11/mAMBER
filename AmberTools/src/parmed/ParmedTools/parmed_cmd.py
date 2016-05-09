"""
This sets up the command interpreter for textual ParmEd (parmed.py).
"""

# Load some system modules that may be useful for various users in shell mode
import math
from chemistry.amber.readparm import AmberParm

import cmd
from ParmedTools import ParmedActions
from ParmedTools.argumentlist import ArgumentList
from ParmedTools.exceptions import InterpreterError

class ParmedCmd(cmd.Cmd):
   """
   ParmEd command interpreter. The docstrings for each do_* command are simply
   copied from the docstring for that Action's docstring
   """

   prompt = "> "

   def __init__(self, parm, stdin=None, stdout=None):
      """ Load a topology file into the interpreter """
      self.parm = parm
      cmd.Cmd.__init__(self, stdin=stdin, stdout=stdout)

   def emptyline(self):
      """ Override emptyline so that empty lines are simply ignored """
      pass

   def precmd(self, line):
      """ Override this to strip out comments, but preserve #s in quotes """
      in_quotes = False
      quote_type = ''
      idx = -1
      for i, char in enumerate(line):
         if in_quotes:
            if char == quote_type:
               in_quotes = False
               continue
            continue
         if char == '"':
            in_quotes = True
            quote_type = '"'
            continue
         if char == "'":
            in_quotes = True
            quote_type = "'"
            continue
         if char == '#':
            idx = i
            break

      if idx < 0:
         return line
      else:
         return line[:idx]

   def parseline(self, line):
      """
      Override parseline so that it will set args as ArgumentList
      """
      line = line.strip()
      if not line:
         return None, None, line
      elif line[0] == '?':
         line = 'help ' + line[1:]
      elif line[0] == '!':
         if hasattr(self, 'do_shell'):
            line = 'shell ' + line[1:]
         else:
            return None, None, line
      i, n = 0, len(line)
      while i < n and line[i] in self.identchars:
         i = i+1
      cmd, arg = line[:i], ArgumentList(line[i:].strip())
      if len(line) == 4 and line == 'help':
         arg = ''
      if line[:5] == 'help ':
         arg = line[i:].strip()
      return cmd, arg, line

   def do_parmout(self, line):
      # Store this action for later use.  This action is unique in that respect
      self.parmout = ParmedActions.parmout(self.parm, line)
   
   def _normaldo(self, ActionClass, line):
      """ The standard action command does this stuff """
      action = ActionClass(self.parm, line)
      self.stdout.write('%s\n' % action)
      action.execute()

   def do_setOverwrite(self, line):
      self._normaldo(ParmedActions.setoverwrite, line)

   def do_writeFrcmod(self, line):
      self._normaldo(ParmedActions.writefrcmod, line)

   def do_loadRestrt(self, line):
      self._normaldo(ParmedActions.loadrestrt, line)

   def do_writeOFF(self, line):
      self._normaldo(ParmedActions.writeoff, line)

   def do_changeRadii(self, line):
      self._normaldo(ParmedActions.changeradii, line)

   def do_changeLJPair(self, line):
      self._normaldo(ParmedActions.changeljpair, line)

   def do_changeLJ14Pair(self, line):
      self._normaldo(ParmedActions.changelj14pair, line)

   def do_checkValidity(self, line):
      self._normaldo(ParmedActions.checkvalidity, line)

   def do_change(self, line):
      self._normaldo(ParmedActions.change, line)

   def do_printInfo(self, line):
      self._normaldo(ParmedActions.printinfo, line)

   def do_addLJType(self, line):
      self._normaldo(ParmedActions.addljtype, line)

   def do_outparm(self, line):
      self._normaldo(ParmedActions.outparm, line)

   def do_printLJTypes(self, line):
      self._normaldo(ParmedActions.printljtypes, line)

   def do_scee(self, line):
      self._normaldo(ParmedActions.scee, line)

   def do_scnb(self, line):
      self._normaldo(ParmedActions.scnb, line)

   def do_changeLJSingleType(self, line):
      self._normaldo(ParmedActions.changeljsingletype, line)

   def do_printDetails(self, line):
      self._normaldo(ParmedActions.printdetails, line)

   def do_printFlags(self, line):
      self._normaldo(ParmedActions.printflags, line)

   def do_printPointers(self, line):
      self._normaldo(ParmedActions.printpointers, line)

   def do_setMolecules(self, line):
      self._normaldo(ParmedActions.setmolecules, line)

   def do_combineMolecules(self, line):
      self._normaldo(ParmedActions.combinemolecules, line)

#   def do_addCoarseGrain(self, line):
#      self._normaldo(ParmedActions.addcoarsegrain, line)

   def do_changeProtState(self, line):
      self._normaldo(ParmedActions.changeprotstate, line)

   def do_netCharge(self, line):
      self._normaldo(ParmedActions.netcharge, line)

   def do_strip(self, line):
      self._normaldo(ParmedActions.strip, line)

   def do_defineSolvent(self, line):
      self._normaldo(ParmedActions.definesolvent, line)

   def do_addExclusions(self, line):
      self._normaldo(ParmedActions.addexclusions, line)

   def do_printBonds(self, line):
      self._normaldo(ParmedActions.printbonds, line)

   def do_printAngles(self, line):
      self._normaldo(ParmedActions.printangles, line)

   def do_printDihedrals(self, line):
      self._normaldo(ParmedActions.printdihedrals, line)

   def do_setBond(self, line):
      self._normaldo(ParmedActions.setbond, line)

   def do_setAngle(self, line):
      self._normaldo(ParmedActions.setangle, line)

   def do_addDihedral(self, line):
      self._normaldo(ParmedActions.adddihedral, line)

   def do_addAtomicNumber(self, line):
      self._normaldo(ParmedActions.addatomicnumber, line)

   def do_deleteDihedral(self, line):
      self._normaldo(ParmedActions.deletedihedral, line)

   def do_printLJMatrix(self, line):
      self._normaldo(ParmedActions.printljmatrix, line)

   def do_tiMerge(self, line):
      self._normaldo(ParmedActions.timerge, line)

   def do_source(self, line):
      action = ParmedActions.source(self.parm, line)
      self.stdout.write('%s\n' % action)
      _cmd = ParmedCmd(self.parm, stdin=open(action.filename, 'r'),
                       stdout=self.stdout)
      _cmd.prompt = ''
      _cmd.interpreter = self.interpreter
      _cmd.use_rawinput = 0
      _cmd.cmdloop()
      
   def do_listParms(self, line):
      self._normaldo(ParmedActions.listparms, line)

   def do_ls(self, line):
      self._normaldo(ParmedActions.ls, line)

   def do_parm(self, line):
      self._normaldo(ParmedActions.parm, line)

   def do_go(self, line):
      """
      Stops reading commands and executes any 'parmout' command that had been
      issued
      """
      if hasattr(self, 'parmout'):
         self.stdout.write('%s\n' % self.parmout)
         self.parmout.execute()

      return True

   # EOF is treated the same as "go"
   do_EOF = do_go

   def do_quit(self, line):
      """
      Stops reading commands and quits WITHOUT running the last 'parmout' 
      command that had been issued
      """
      return True

   def default(self, line):
      """ Default behavior """
      mycmd = line.split()[0].lower()
      if mycmd in ('lawsuit', 'math', 'warnings') or not \
            hasattr(ParmedActions, mycmd):
         self.stdout.write("%s command not recognized\n" % line.split()[0])
      if hasattr(ParmedActions, mycmd):
         args = ArgumentList(line.split()[1:])
         self._normaldo(getattr(ParmedActions, mycmd), args)

   def do_shell(self, line):
      """ Support the limited interpreter """
      if not self.interpreter:
         raise InterpreterError("Interpreter not enabled! Use '-e' to enable")
      line = str(line)
      amber_prmtop = self.parm
      if line.strip() == '!':
         self._python_shell()
      else:
         try:
            exec(line.strip())
         except BaseException, err:
            self.stdout.write("%s: %s\n" % (type(err).__name__, err))

   def _python_shell(self):
      """ Drop into a limited interpreter and read until we see !! """
      from ParmedTools import ParmedActions
      python_interpreter = PythonCmd(stdin=self.stdin, stdout=self.stdout)
      python_interpreter.use_rawinput = self.use_rawinput
      if not self.prompt: python_interpreter.prompt = ''
      python_interpreter.cmdloop()
      try:
         amber_prmtop = self.parm
         exec(python_interpreter.command_string)
      except BaseException, err:
         self.stdout.write("%s: %s\n" % (type(err).__name__, err))
      
   def do_help(self, arg):
      """ Modify the original do_help to pull docstrings from ParmedActions """
      if arg:
         # XXX check arg syntax
         try:
            func = getattr(self, 'help_' + arg)
         except AttributeError:
            try:
               doc=getattr(self, 'do_' + arg).__doc__
               if doc:
                  self.stdout.write("%s\n" % str(doc))
                  return
            except AttributeError:
               pass
            try:
               _action = getattr(ParmedActions, arg.lower())
               if _action.needs_parm:
                  doc = '%s [parm <idx>|<name>]\n%s'
               else:
                  doc = '%s\n%s'
               doc = doc % (ParmedActions.Usages[arg.lower()],
                            getattr(ParmedActions, arg.lower()).__doc__)
               if doc:
                  self.stdout.write("%s\n"%str(doc))
                  return
            except (AttributeError, KeyError):
               pass
            self.stdout.write("%s\n" % str(self.nohelp % (arg,)))
            return
            func()
      else:
         names = self.get_names()
         cmds_doc = []
         cmds_undoc = []
         help = {}
         for name in names:
            if name[:5] == 'help_':
               help[name[5:]]=1
         names.sort()
         # There can be duplicates if routines overridden
         prevname = ''
         for name in names:
            if name[:3] == 'do_':
               if name == prevname:
                  continue
               prevname = name
               cmd=name[3:]
               if cmd in help:
                  cmds_doc.append(cmd)
                  del help[cmd]
               elif getattr(self, name).__doc__:
                  cmds_doc.append(cmd)
               elif hasattr(ParmedActions, cmd.lower()) and \
                     cmd.lower() in ParmedActions.Usages.keys():
                  cmds_doc.append(cmd)
               else:
                  cmds_undoc.append(cmd)
         self.stdout.write("%s\n"%str(self.doc_leader))
         self.print_topics(self.doc_header,   cmds_doc,   15,80)
         self.print_topics(self.misc_header,  help.keys(),15,80)
         self.print_topics(self.undoc_header, cmds_undoc, 15,80)

class PythonCmd(cmd.Cmd):
   """ Command interpreter for limited Python interpreter """
   prompt = "py >>> "

   def __init__(self, stdin=None, stdout=None):
      cmd.Cmd.__init__(self, stdin=stdin, stdout=stdout)
      self.command_string = ""

   def emptyline(self):
      """ Ignore all empty lines """
      self.command_string += "\n"

   def do_shell(self, line):
      """ Break out of the shell """
      return True

   def do_EOF(self, line):
      """ Break out of the shell """
      raise InterpreterError("EOF hit while in interpreter!")

   def parseline(self, line):
      """
      Override parseline
      """
      if not line.strip():
         return False, False, False
      return line, line, line

   def default(self, line):
      """ Add this onto the command string """
      if line.strip() == '!!': return True
      self.command_string += line + '\n'
