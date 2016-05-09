#! PYTHONEXE

# Load system modules.
import signal
import sys
import os
from argparse import ArgumentParser

# Load custom modules
from ParmedTools.logos import Logo
from ParmedTools.parmed_cmd import ParmedCmd
from ParmedTools.parmlist import ParmList
from ParmedTools import __version__

# Set up new excepthook to clean up fatal exception printouts
debug = False

def interrupted(*args, **kwargs):
   """ Handle interruptions gracefully """
   sys.stdout.write('Interrupted\n')
   sys.exit(1)

signal.signal(signal.SIGINT, interrupted)

def excepthook(exception_type, exception_value, tb):
   """ Replaces excepthook so we can control the printing of 
       tracebacks. Those are helpful for debugging purposes, but 
       may be unsightly to users. debug set above controls this
       behavior, and a CL flag allows users to set this.
   """
   import traceback
   if debug: traceback.print_tb(tb)
   sys.stderr.write('%s: %s\n' % (exception_type.__name__, exception_value))

sys.excepthook = excepthook

# Set up parser
parser = ArgumentParser(version="%%(prog)s: Version %s" % __version__)
group = parser.add_argument_group('Input Files')
group.add_argument('-i', '--input', dest='script', default=[], 
                   metavar='FILE', help='''Script with ParmEd commands to
                   execute. Default reads from stdin. Can be specified multiple
                   times to process multiple input files.''', action='append')
group.add_argument('-p', '--parm', dest='prmtop', default=[],
                   metavar='<prmtop>', action='append', help='''List of topology
                   files to load into ParmEd. Can be specified multiple times to
                   process multiple topologies.''')
group = parser.add_argument_group('ParmEd Operating Options')
parser.add_argument('-d', '--debug', dest='debug', 
                    action='store_true', default=False,
                    help='''Show tracebacks for any uncaught exceptions. Useful
                    for debugging. OFF by default''')
parser.add_argument('-e', '--enable-interpreter', dest='interpreter',
                    action='store_true', default=False,
                    help='''Print how each action is parsed and show which line
                    of which file an error occurred on (verbose tracebacks).
                    OFF by default.''')
parser.add_argument('-q', '--quiet', dest='debug', action='store_false',
                    help='Disable verbose tracebacks. Reverses -d/--debug')
parser.add_argument('--prompt', dest='prompt', default='>', metavar='PROMPT',
                    help='String to use as a command prompt.')
parser.add_argument('-n', '--no-splash', dest='printlogo', action='store_false',
                    help='Prevent printing the greeting logo.', default=True)
parser.add_argument('prmtop_cl', nargs='?', metavar='<prmtop>', default=None,
                    help='Topology file to analyze.')
parser.add_argument('script_cl', nargs='?', metavar='<script>', default=None,
                    help='Script with ParmEd commands to execute.')

opt = parser.parse_args()

debug = opt.debug

# If the user specified a prmtop and script in the 'old' way, append them to the
# list constructed via the --parm and -i flags -- they come at the end
if opt.script_cl is not None:
   opt.script.append(opt.script_cl)

if opt.prmtop_cl is not None:
   opt.prmtop.append(opt.prmtop_cl)

# Load the splash screen
if opt.printlogo:
   splash = Logo()
   print splash

amber_prmtop = ParmList()
for parm in opt.prmtop:
   amber_prmtop.add_parm(parm)
   print 'Loaded Amber topology file %s\n' % parm

try:
   # Clear this out of the namespace for scripting purposes (if present)
   del parm
except NameError:
   pass

if len(opt.script) > 0:
   # Read from the list of scripts
   print opt.script
   for script in opt.script:
      if not os.path.exists(script):
         print 'Warning: Script file %s cannot be found' % script
      print 'Reading actions from %s\n' % script
      parmed_commands = ParmedCmd(amber_prmtop, stdin=open(script, 'r'))
      parmed_commands.use_rawinput = 0
      parmed_commands.interpreter = opt.interpreter
      parmed_commands.prompt = ''
      # Loop through all of the commands
      parmed_commands.cmdloop()

else:
   parmed_commands = ParmedCmd(amber_prmtop)
   parmed_commands.intro = "Reading input from STDIN..."
   parmed_commands.interpreter = opt.interpreter
   parmed_commands.prompt = opt.prompt.strip() + ' '
   # Loop through all of the commands
   parmed_commands.cmdloop()

print 'Done!'
