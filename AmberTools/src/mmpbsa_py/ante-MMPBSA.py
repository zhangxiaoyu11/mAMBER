#! PYTHONEXE

from MMPBSA_mods.findprogs import which
from optparse import OptionParser
from subprocess import Popen, PIPE
import sys

clparser = OptionParser()
clparser.add_option('-p', '--prmtop', dest='prmtop', default=None,
         help='Input "dry" complex topology or solvated complex topology')
clparser.add_option('-c', '--complex-prmtop', dest='complex', default=None,
         help='Complex topology file created by stripping PRMTOP of solvent')
clparser.add_option('-r', '--receptor-prmtop', dest='receptor', default=None,
         help='Receptor topology file created by stripping COMPLEX of ligand')
clparser.add_option('-l', '--ligand-prmtop', dest='ligand', default=None,
         help='Ligand topology file created by stripping COMPLEX of receptor')
clparser.add_option('-s', '--strip-mask', dest='strip_mask', default=None,
         help='Amber mask of atoms needed to be stripped from PRMTOP to make ' +
              'the COMPLEX topology file')
clparser.add_option('-m', '--receptor-mask', dest='receptor_mask', default=None,
         help='Amber mask of atoms needed to be stripped from COMPLEX to ' +
              'create RECEPTOR. Cannot specify with -n/--ligand-mask')
clparser.add_option('-n', '--ligand-mask', dest='ligand_mask', default=None,
         help='Amber mask of atoms needed to be stripped from COMPLEX to ' +
              'create LIGAND. Cannot specify with -m/--receptor-mask')
clparser.add_option('--radii', dest='radius_set', default=None,
         help='PB/GB Radius set to set in the generated topology files. ' +
              'This is equivalent to "set PBRadii <radius>" in LEaP. ' +
              'Options are bondi, mbondi2, mbondi3, amber6, and mbondi ' +
              'and the default is to use the existing radii.')
(opt, args) = clparser.parse_args()

# Check for illegal CL options
if not opt.prmtop:
   clparser.print_help()
   sys.exit(1)

if len(args):
   print >> sys.stderr, 'Error: Unknown command-line arguments: %s' % (
            ', '.join(args))
   sys.exit(1)

if opt.receptor_mask and opt.ligand_mask:
   print >> sys.stderr, 'Command-line Error: Cannot specify both receptor ' + \
         'and ligand masks!'
   sys.exit(1)

if (opt.receptor_mask or opt.ligand_mask) and not (opt.receptor and opt.ligand):
   print >> sys.stderr, 'Error: You specified ligand or receptor mask, but ' + \
         'not a ligand and receptor topology file!'
   sys.exit(1)

if not (opt.ligand_mask or opt.receptor_mask) and (opt.receptor or opt.ligand):
   print >> sys.stderr, 'Error: You must provide a ligand or receptor mask ' + \
         'for a ligand and receptor topology file!'
   sys.exit(1)

if not (opt.receptor and opt.ligand) and (opt.ligand_mask or opt.receptor_mask):
   print >> sys.stderr, 'Error: You must provide a ligand or receptor mask ' + \
         'for a ligand or receptor topology file.'
   sys.exit(1)

if (opt.receptor and not opt.ligand) or (opt.ligand and not opt.receptor):
   print >> sys.stderr, 'Error: You must specify both ligand and receptor ' + \
         'topologies or neither!'
   sys.exit(1)

if not opt.ligand_mask and not opt.receptor_mask and not opt.strip_mask:
   print >> sys.stderr, 'Error: You did not specify any masks -- I have ' + \
         'nothing to do!'
   sys.exit(1)

if opt.radius_set is not None:
   allowed_radii = ('bondi', 'mbondi', 'mbondi2', 'mbondi3', 'amber6')
   if not opt.radius_set.lower() in allowed_radii:
      print >> sys.stderr, 'Error: Radius set must be one of %s' % (
               ' '.join(allowed_radii))
      sys.exit(1)

# Now load the unspecified mask (just !(specified_mask))
if opt.receptor_mask or opt.ligand_mask:
   if not opt.receptor_mask: opt.receptor_mask = '!(%s)' % opt.ligand_mask
   if not opt.ligand_mask: opt.ligand_mask = '!(%s)' % opt.receptor_mask

parmed = which('parmed.py', search_path=True)

# check that necessary parameters are met
if not parmed:
   print >> sys.stderr, 'Error: parmed.py cannot be found!'
   sys.exit(1)

# Create the stripped topology file
if opt.strip_mask:
   print 'Stripping %s (solvent) from original topology, output is %s' % (
         opt.strip_mask, opt.complex)
   parmed_commands  = "strip %s\n" % opt.strip_mask
   if opt.radius_set is not None:
      parmed_commands += "changeRadii %s\n" % opt.radius_set.lower()
   parmed_commands += "parmout %s\n" % opt.complex
   parmed_commands += "go\n"
   process = Popen([parmed, '-n', '-q', opt.prmtop], stdout=PIPE, stderr=PIPE,
                   stdin=PIPE)
   (output, error) = process.communicate(parmed_commands)
   if process.wait():
      print >> sys.stderr, 'Error: Creating complex topology failed!'
      print >> sys.stderr, output.strip(), '\n', error.strip()
      sys.exit(1)

   print 'Done stripping solvent!\n'
# If we aren't stripping solvent, our complex is our original prmtop
else: opt.complex = opt.prmtop

# Now create the receptor prmtop
if opt.receptor:
   print 'Creating receptor topology file by stripping %s from %s' % (
         opt.ligand_mask, opt.complex)
   
   parmed_commands  = "strip %s\n" % opt.ligand_mask
   if opt.radius_set is not None:
      parmed_commands += "changeRadii %s\n" % opt.radius_set.lower()
   parmed_commands += "parmout %s\n" % opt.receptor
   parmed_commands += "go\n"
   process = Popen([parmed, '-n', '-q', opt.complex], stdout=PIPE, stderr=PIPE,
                   stdin=PIPE)
   (output, error) = process.communicate(parmed_commands)
   if process.wait():
      print >> sys.stderr, 'Error: Creating receptor topology failed!'
      print >> sys.stderr, output.strip(), '\n', error.strip()
      sys.exit(1)

   print 'Done creating receptor topology file!\n'

# Now create the ligand prmtop
if opt.ligand:
   print 'Creating ligand topology file by stripping %s from %s' % (
         opt.receptor_mask, opt.complex)

   parmed_commands  = "strip %s\n" % opt.receptor_mask
   if opt.radius_set is not None:
      parmed_commands += "changeRadii %s\n" % opt.radius_set.lower()
   parmed_commands += "parmout %s\n" % opt.ligand
   parmed_commands += "go\n"
   process = Popen([parmed, '-n', '-q', opt.complex], stdout=PIPE, stderr=PIPE,
                   stdin=PIPE)
   (output, error) = process.communicate(parmed_commands)

   if process.wait():
      print >> sys.stderr, 'Error: Creating ligand topology failed!'
      print >> sys.stderr, output.strip(), '\n', error.strip()
      sys.exit(1)
   
   print 'Done creating ligand topology file!\n'
