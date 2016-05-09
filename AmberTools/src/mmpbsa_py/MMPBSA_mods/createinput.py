"""
This module contains classes for creating input files for various MM/PBSA
calculations.

######################### GPL LICENSE INFO ############################

  Copyright (C) 2009  Dwight McGee, Billy Miller III, and Jason Swails

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

Methods:
   create_inputs -- determines which input files need to be written then 
                    handles it

Classes: 
    NabInput -- base class of mmpbsa_py_energy input files
    GBNabInput -- writes input files for GB (derived from NabInput)
    PBNabInput -- writes input files for PB (derived from NabInput)
    SanderInput -- Base class for sander input files
    SanderGBInput -- writes input files for sander GB (derived from SanderInput)
    SanderPBSAInput- writes input files for sander PB (derived from SanderInput)
    SanderAPBSInput -- writes input files for sander.APBS PB 
                       (derived from SanderInput)
    SanderGBDecomp -- writes GB decomp input files (derived from SanderInput)
    SanderPBDecomp -- writes PB decomp input files (derived from SanderInput)
    QuasiHarmonicInput -- writes a ptraj input file for quasi-harmonic calcs
"""

from MMPBSA_mods.exceptions import CreateInputError
from copy import deepcopy

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def create_inputs(INPUT, prmtop_system, pre):
   """ Creates the input files for all necessary calculations """
   stability = prmtop_system.stability

   # First check if we are running decomp
   if INPUT['decomprun']:
      # Get the cards that will go in the complex mdin file
      com_card, rc_card, lc_card = prmtop_system.Group(INPUT['print_res'],True)
      junk, rec_card, lig_card = prmtop_system.Group(INPUT['print_res'],False,
                                 rec_group_type='RES', lig_group_type='RES')
      full_com, full_rc, full_lc = prmtop_system.Group('all',True)
      junk, full_rec, full_lig = prmtop_system.Group('all',False)
      
      # Convert the strings into card objects
      # Now create the mdin objects
      if INPUT['gbrun']:
         if stability:
            rec_res = ['Residues considered as REC', full_rc]
            pri_res = ['Residues to print', com_card]
            com_mdin = SanderGBDecomp(INPUT, rec_res, pri_res)
            com_mdin.write_input(pre + 'gb_decomp_com.mdin')
         else:
            rec_res = ['Residues considered as REC', full_rc]
            lig_res = ['Residues considered as LIG', full_lc]
            pri_res = ['Residues to print', com_card]
            com_mdin = SanderGBDecomp(INPUT, rec_res, lig_res, pri_res)
            rec_res = ['Residues considered as REC', full_rec]
            pri_res = ['Residues to print', rec_card]
            rec_mdin = SanderGBDecomp(INPUT, rec_res, pri_res)
            lig_res = ['Residues considered as LIG', full_lig]
            pri_res = ['Residues to print', lig_card]
            lig_mdin = SanderGBDecomp(INPUT, lig_res, pri_res)
            com_mdin.write_input(pre + 'gb_decomp_com.mdin')
            rec_mdin.write_input(pre + 'gb_decomp_rec.mdin')
            lig_mdin.write_input(pre + 'gb_decomp_lig.mdin')
      if INPUT['pbrun']:
         if stability:
            rec_res = ['Residues considered as REC', full_rc]
            pri_res = ['Residues to print', com_card]
            com_mdin = SanderPBDecomp(INPUT, rec_res, pri_res)
            com_mdin.write_input(pre + 'pb_decomp_com.mdin')
         else:
            rec_res = ['Residues considered as REC', full_rc]
            lig_res = ['Residues considered as LIG', full_lc]
            pri_res = ['Residues to print', com_card]
            com_mdin = SanderPBDecomp(INPUT, rec_res, lig_res, pri_res)
            rec_res = ['Residues considered as REC', full_rec]
            pri_res = ['Residues to print', rec_card]
            rec_mdin = SanderPBDecomp(INPUT, rec_res, pri_res)
            lig_res = ['Residues considered as LIG', full_lig]
            pri_res = ['Residues to print', lig_card]
            lig_mdin = SanderPBDecomp(INPUT, lig_res, pri_res)
            com_mdin.write_input(pre + 'pb_decomp_com.mdin')
            rec_mdin.write_input(pre + 'pb_decomp_rec.mdin')
            lig_mdin.write_input(pre + 'pb_decomp_lig.mdin')

   else: # not decomp
      
      if INPUT['gbrun']:
         gb_prog = 'mmpbsa_py_energy'
         if INPUT['use_sander'] or INPUT['igb'] == 8 or INPUT['ifqnt'] == 1:
            gb_prog = 'sander'

         if gb_prog == 'mmpbsa_py_energy':
            gb_mdin = GBNabInput(INPUT)
         else:
            # We need separate input files for QM/MMGBSA
            if INPUT['ifqnt']:
               com_input = deepcopy(INPUT)
               rec_input = deepcopy(INPUT)
               lig_input = deepcopy(INPUT)
               (com_input['qmmask'], rec_input['qmmask'], 
                lig_input['qmmask']) = prmtop_system.Mask(INPUT['qm_residues'],
                                        in_complex=False)
               if not com_input['qmmask']:
                  raise CreateInputError('No valid QM residues chosen!')
               com_input['qm_theory'] = "'%s'" % com_input['qm_theory']
               com_input['qmmask'] = "'%s'" % com_input['qmmask']
               gb_mdin = SanderGBInput(com_input)
               gb_mdin.write_input(pre + 'gb_qmmm_com.mdin')
               if not stability:
                  if not rec_input['qmmask']: rec_input['ifqnt'] = 0
                  else: rec_input['qmmask'] = "'%s'" % rec_input['qmmask']
                  rec_input['qm_theory'] = "'%s'" % rec_input['qm_theory']
                  gb_mdin = SanderGBInput(rec_input)
                  gb_mdin.write_input(pre + 'gb_qmmm_rec.mdin')
                  if not lig_input['qmmask']: lig_input['ifqnt'] = 0
                  else: lig_input['qmmask'] = "'%s'" % lig_input['qmmask']
                  lig_input['qm_theory'] = "'%s'" % lig_input['qm_theory']
                  gb_mdin = SanderGBInput(lig_input)
                  gb_mdin.write_input(pre + 'gb_qmmm_lig.mdin')
            else:
               gb_mdin = SanderGBInput(INPUT)

         gb_mdin.write_input(pre + 'gb.mdin')

      if INPUT['pbrun']:
         pb_prog = 'mmpbsa_py_energy'
         if INPUT['sander_apbs']:
            pb_prog = 'sander.APBS'
         elif INPUT['use_sander']:
            pb_prog = 'sander'

         if pb_prog == 'mmpbsa_py_energy':
            pb_mdin = PBNabInput(INPUT)
         elif pb_prog == 'sander.APBS':
            pb_mdin = SanderAPBSInput(INPUT)
         else:
            pb_mdin = SanderPBSAInput(INPUT)
         
         pb_mdin.write_input(pre + 'pb.mdin')
         
   #end if decomprun

   if INPUT['entropy']: # quasi-harmonic approximation input file
      trj_suffix = 'mdcrd'
      if INPUT['netcdf']: trj_suffix = 'nc'
      com_mask, rec_mask, lig_mask = prmtop_system.Mask('all', True)
      if not INPUT['mutant_only']:
         qh_in = QuasiHarmonicInput(com_mask, rec_mask, lig_mask,
                  stability=stability, prefix=pre, trj_suffix=trj_suffix)
         qh_in.write_input(pre + 'ptrajentropy.in')
      if INPUT['alarun']:
         qh_in = QuasiHarmonicInput(com_mask, rec_mask, lig_mask,
                  stability=stability, prefix=pre + 'mutant_',
                  trj_suffix=trj_suffix)
         qh_in.write_input(pre + 'mutant_ptrajentropy.in')

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class NabInput(object):
   """ 
   Creates an input file for the Nab energy program. It must follow the format
   expected by the nab program:

   GB (or PB)
   variable1 = value1
   variable2 = value2
   ...
   variableN = valueN

   input_items is a dictionary of recognized keywords in the input file

   name_map is a dictionary that maps each key from input_items to the 
   corresponding key from the INPUT dictionary in MMPBSA.py
   """
   input_items = {'foo' : 'bar'}  # This should be overridden in derived classes
   name_map    = {'foo' : 'orig'} # This should be overridden in derived classes
   calculation_type = 'GB'        # This should be overridden in derived classes

   def __init__(self, INPUT):
      # Replace each of my input_items with the value held in INPUT if applicable
      # but if it's not in INPUT, don't pitch a fit (also ignore the possibility
      # that not every key is mapped in name_map, since only those that are in
      # INPUT will be mapped
      for key in self.input_items.keys():
         try: self.input_items[key] = INPUT[self.name_map[key]]
         except KeyError: pass

   def write_input(self, filename):
      """ Writes the input file to filename """
      if isinstance(filename, str):
         infile = open(filename, 'w')
      elif isinstance(filename, file):
         infile = filename
      else:
         raise TypeError('write_input expects a str file name or file object')

      infile.write('%s\n' % self.calculation_type)
      for key in self.input_items.keys():
         infile.write('%s = %s\n' % (key, self.input_items[key]))
      infile.close()

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class GBNabInput(NabInput):
   """ The GB input file for the mmpbsa_py_energy NAB program """
   input_items = {'igb'     : 1,
                  'extdiel' : 78.3,
                  'saltcon' : 0.0,
                  'surften' : 0.0072,
                  'rgbmax'  : 25.0 }
   name_map = {'igb' : 'igb', 
               'extdiel' : 'extdiel',
               'saltcon' : 'saltcon',
               'surften' : 'surften'}
   calculation_type = 'GB'

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class PBNabInput(NabInput):
   """ The PB input file for the mmpbsa_py_energy NAB program """
   input_items = {'inp'       : 2, 'smoothopt' : 1, 'radiopt'   : 1,
                  'npbopt'    : 0, 'solvopt'   : 1, 'maxitn'    : 1000,
                  'nfocus'    : 2, 'bcopt'     : 5, 'eneopt'    : 2,
                  'fscale'    : 8, 'epsin'     : 1.0, 'epsout'    : 80.0,
                  'istrng'    : 0.0, 'dprob'   : 1.4, 'iprob'     : 2.0,
                  'accept'    : 0.001, 'fillratio' : 4, 'space'     : 0.5,
                  'cutnb'     : 0, 'sprob'     : 0.557, 'dbfopt'    : 1,
                  'cavity_surften' : 0.0378, 'cavity_offset'  : -0.5692 }
   name_map = {'inp' : 'inp', 'radiopt' : 'radiopt', 'maxitn' : 'linit',
               'epsin' : 'indi', 'epsout' : 'exdi', 'istrng' : 'istrng',
               'dprob' : 'prbrad', 'fillratio' : 'fillratio', 'space' : 'scale',
               'cavity_surften' : 'cavity_surften',
               'cavity_offset' : 'cavity_offset' }
   calculation_type = 'PB'

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class SanderInput(object):
   """ Base class sander input file """
   program = 'sander'                         # This runs with sander
   input_items = {'foo' : 'bar'}              # replace this in derived classes
   name_map = {'foo' : 'orig'}                # replace this in derived classes
   parent_namelist = {'foo' : 'foo_namelist'} # replace this in derived classes
   
   def __init__(self, INPUT):
      from mdin import mdin
      self.mdin = mdin(self.program)
      self.mdin.title = 'File generated by MMPBSA.py'
      for key in self.input_items.keys():
         try: self.mdin.change(self.parent_namelist[key], key,
                               INPUT[self.name_map[key]])
         except KeyError:
            self.mdin.change(self.parent_namelist[key], key,
                             self.input_items[key])
      self.mdin.change('cntrl', 'ioutfm', int(bool(INPUT['netcdf'])))

   def write_input(self, filename):
      """ Write the mdin file """
      if not isinstance(filename, str):
         raise TypeError('Sander input requires a str filename')
      self.mdin.write(filename)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class SanderGBInput(SanderInput):
   """ GB sander input file """
   input_items = {'ntb':0,    'cut':999.0,    'nsnb':99999,
                  'imin':5,   'maxcyc':1,     'ncyc':0,
                  'igb':5,    'saltcon':0.0,  'intdiel':1.0,
                  'gbsa':0,   'extdiel':80.0, 'surften':0.0072,
                  'ioutfm':0, 'idecomp':0,    'offset':-999999.0,
                  'dec_verbose':0, 'ifqnt':0, 'qmmask':'',
                  'qm_theory':'', 'qmcharge':0, 'qmgb':2,
                  'qmcut':999.0 }
   parent_namelist = {'ntb':'cntrl',     'cut':'cntrl',     'nsnb':'cntrl',
                      'imin':'cntrl',    'maxcyc':'cntrl',  'ncyc':'cntrl',
                      'igb':'cntrl',     'saltcon':'cntrl', 'intdiel':'cntrl',
                      'extdiel':'cntrl', 'gbsa':'cntrl',    'surften':'cntrl',
                      'offset':'cntrl',  'idecomp':'cntrl', 'ioutfm':'cntrl',
                      'dec_verbose':'cntrl', 'ifqnt':'cntrl', 'qmmask':'qmmm',
                      'qm_theory':'qmmm', 'qmcharge':'qmmm', 'qmgb':'qmmm',
                      'qmcut':'qmmm'}
   name_map = {'igb':'igb',  'saltcon':'saltcon', 'intdiel':'intdiel',
               'extdiel':'extdiel', 'gbsa':'gbsa', 'surften':'surften',
               'idecomp':'idecomp', 'ioutfm':'netcdf', 'qmmask':'qmmask',
               'qm_theory':'qm_theory', 'ifqnt':'ifqnt', 'qmcharge':'qmcharge',
               'qmcut':'qmcut', 'dec_verbose':'dec_verbose'}

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class SanderPBSAInput(SanderInput):
   """ PB sander input file """
   input_items = {'ntb':0, 'cut':999.0, 'nsnb':99999,
                  'imin':5, 'maxcyc':1, 'igb':10,
                  'ipb':2, 'inp':2, 'idecomp':0,
                  'dbfopt':1, 'epsin':1.0, 'epsout':80.0,
                  'istrng':0.0, 'radiopt':1, 'sprob':0.557, 'dprob':1.4,
                  'space':0.5, 'maxitn':1000, 'cavity_surften':0.0378,
                  'cavity_offset':-0.5692, 'fillratio':4.0,
                  'ioutfm':0, 'dec_verbose':0}
   name_map = {'inp':'inp', 'idecomp':'idecomp', 'epsin':'indi',
               'epsout':'exdi', 'istrng':'istrng', 'radiopt':'radiopt',
               'dprob':'prbrad', 'space':'scale', 'maxitn':'linit',
               'cavity_surften':'cavity_surften', 
               'cavity_offset':'cavity_offset', 'fillratio':'fillratio',
               'ioutfm':'netcdf', 'dec_verbose':'dec_verbose'}
   parent_namelist = {'ntb':'cntrl', 'cut':'cntrl', 'nsnb':'cntrl',
                      'imin':'cntrl', 'maxcyc':'cntrl', 'igb':'cntrl',
                      'ipb':'cntrl', 'inp':'cntrl', 'idecomp':'cntrl',
                      'dbfopt':'pb', 'epsin':'pb', 'epsout':'pb',
                      'istrng':'pb', 'radiopt':'pb', 'sprob':'pb', 'dprob':'pb',
                      'space':'pb', 'maxitn':'pb', 'cavity_offset':'pb',
                      'cavity_surften':'pb', 'fillratio':'pb', 
                      'ioutfm':'cntrl', 'dec_verbose':'cntrl'}
   def __init__(self, INPUT):
      # We need to change istrng to mM (from M).
      SanderInput.__init__(self, INPUT)
      self.mdin.change('pb', 'istrng', INPUT['istrng']*1000)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class SanderAPBSInput(SanderInput):
   """ PB sander input file using APBS as the PB solver """
   program = 'sander.APBS'
   input_items = {'ntb':0, 'cut':999.0, 'nsnb':99999,
                  'imin':5, 'maxcyc':1, 'igb':6,
                  'apbs_print':1, 'calc_type':0, 'cmeth':1,
                  'bcfl':2, 'srfm':2, 'chgm':1, 'pdie':1.0,
                  'sdie':80.0, 'srad':1.4, 'nion':2, 'ionq':'1.0,-1.0',
                  'ionc':'0.0,0.0', 'ionrr':'2.0,2.0', 'radiopt':0,
                  'calcforce':0, 'calcnpenergy':1, 'grid':'0.5,0.5,0.5',
                  'gamma':0.00542, 'ioutfm':0}
   name_map = {'inp':'inp', 'idecomp':'idecomp', 'pdie':'indi',
               'sdie':'exdi', 'ionc':'istrng', 'radiopt':'radiopt',
               'srad':'prbrad', 'grid':'scale', 'gamma':'cavity_surften',
               'ioutfm':'netcdf'}
   parent_namelist = {'ntb':'cntrl', 'cut':'cntrl', 'nsnb':'cntrl',
                      'imin':'cntrl', 'maxcyc':'cntrl', 'igb':'cntrl',
                      'ipb':'cntrl', 'inp':'cntrl', 'idecomp':'cntrl',
                      'ioutfm':'cntrl', 'apbs_print':'apbs', 'calc_type':'apbs',
                      'cmeth':'apbs', 'bcfl':'apbs', 'srfm':'apbs', 
                      'chgm':'apbs', 'pdie':'apbs', 'sdie':'apbs', 
                      'srad':'apbs', 'nion':'apbs', 'ionq':'apbs', 
                      'ionc':'apbs', 'ionrr':'apbs', 'radiopt':'apbs', 
                      'calcforce':'apbs', 'grid':'apbs', 'gamma':'apbs',
                      'calcnpenergy':'apbs'}
   def __init__(self, INPUT):
      SanderInput.__init__(self, INPUT)
      # We also have to make some modifications specific to sander.APBS
      self.mdin.change('apbs', 'grid', '%s,%s,%s' % (INPUT['scale'],
         INPUT['scale'], INPUT['scale']))
      self.mdin.change('apbs', 'ionc', '%s,%s' % (INPUT['istrng'], 
         INPUT['istrng']))
      self.mdin.change('apbs', 'gamma', INPUT['cavity_surften']*4.184)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class SanderGBDecomp(SanderGBInput):
   """ GB decomposition input file for sander. In addition to the INPUT dictionary,
       this class also needs several GROUP cards to define the 
   """
   def __init__(self, INPUT, *cards):
      SanderGBInput.__init__(self, INPUT)
      for card in cards:
         self.mdin.AddCard(card[0], card[1])

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class SanderPBDecomp(SanderPBSAInput):
   """ PB decomposition input file for sander """
   def __init__(self, INPUT, *cards):
      SanderPBSAInput.__init__(self, INPUT)
      for card in cards:
         self.mdin.AddCard(card[0], card[1])

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class QuasiHarmonicInput(object):
   """ Write a ptraj input file to do a quasi-harmonic entropy calculation """

   def __init__(self, com_mask, rec_mask, lig_mask, stability=False,
                prefix='_MMPBSA_', trj_suffix='mdcrd'):
      self.file_string = ("trajin %scomplex.%s\n" % (prefix, trj_suffix) +
                          "reference %savgcomplex.pdb\n" % (prefix) +
                          "rms mass reference %s\n" % (com_mask) +
                          "matrix mwcovar name comp.matrix %s\n" % (com_mask) +
                          "analyze matrix comp.matrix out " +
                          "%scomplex_entropy.out thermo reduce\n" % (prefix))
      if not stability:
         self.file_string = (self.file_string +
                           "rms mass reference %s\n" % (rec_mask) +
                           "matrix mwcovar name rec.matrix %s\n" % (rec_mask) +
                           "analyze matrix rec.matrix out " +
                           "%sreceptor_entropy.out thermo reduce\n" % (prefix) +
                           "rms mass reference %s\n" % (lig_mask) +
                           "matrix mwcovar name lig.matrix %s\n" % (lig_mask) +
                           "analyze matrix lig.matrix out " +
                           "%sligand_entropy.out thermo reduce\n" % (prefix))

   def write_input(self, filename):
      """ Writes the input file """
      if not isinstance(filename, str):
         raise TypeError('QuasiHarmonicInput expects str type for file name')
      infile = open(filename, 'w')
      infile.write(self.file_string)
      infile.close()

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
