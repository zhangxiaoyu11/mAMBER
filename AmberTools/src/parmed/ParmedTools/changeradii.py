from ParmedTools.exceptions import ChangeRadiiError

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def bondi(parm):
   """ Sets the bondi radii """
   molecule = parm.ToMolecule()
   for i in range(parm.pointers['NATOM']):
      # Radius of C atom depends on what type it is
      if molecule.elements[i] == 'C':
         if parm.parm_data['AMBER_ATOM_TYPE'][i][:2] in ['C1', 'C2', 'C3']:
            parm.parm_data['RADII'][i] = 2.2
         else:
            parm.parm_data['RADII'][i] = 1.7
      # All other elements have fixed radii for all types/partners
      elif molecule.elements[i] == 'H': parm.parm_data['RADII'][i] = 1.2
      elif molecule.elements[i] == 'N': parm.parm_data['RADII'][i] = 1.55
      elif molecule.elements[i] == 'O': parm.parm_data['RADII'][i] = 1.5
      elif molecule.elements[i] == 'F': parm.parm_data['RADII'][i] = 1.5
      elif molecule.elements[i] == 'Si': parm.parm_data['RADII'][i] = 2.1
      elif molecule.elements[i] == 'P': parm.parm_data['RADII'][i] = 1.85
      elif molecule.elements[i] == 'S': parm.parm_data['RADII'][i] = 1.8
      elif molecule.elements[i] == 'Cl': parm.parm_data['RADII'][i] = 1.5
      else: parm.parm_data['RADII'][i] = 1.5

   parm.parm_data['RADIUS_SET'][0] = 'Bondi radii (bondi)'
   _screen1(parm, molecule)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def amber6(parm):
   """ Sets the amber6 radii """
   molecule = parm.ToMolecule()
   for i in range(parm.pointers['NATOM']):
      # Radius of H atom depends on element it is bonded to
      if molecule.elements[i] == 'H':
         if molecule.elements[molecule.bonds[i][0]] in ['C']:
            parm.parm_data['RADII'][i] = 1.3
         elif molecule.elements[molecule.bonds[i][0]] in ['O', 'S']:
            parm.parm_data['RADII'][i] = 0.8
         else:
            parm.parm_data['RADII'][i] = 1.2
      # Radius of C atom depends on what type it is
      elif molecule.elements[i] == 'C':
         if parm.parm_data['AMBER_ATOM_TYPE'][i][:2] in ['C1', 'C2', 'C3']:
            parm.parm_data['RADII'][i] = 2.2
         else:
            parm.parm_data['RADII'][i] = 1.7
      # All other elements have fixed radii for all types/partners
      elif molecule.elements[i] == 'N': parm.parm_data['RADII'][i] = 1.55
      elif molecule.elements[i] == 'O': parm.parm_data['RADII'][i] = 1.5
      elif molecule.elements[i] == 'F': parm.parm_data['RADII'][i] = 1.5
      elif molecule.elements[i] == 'Si': parm.parm_data['RADII'][i] = 2.1
      elif molecule.elements[i] == 'P': parm.parm_data['RADII'][i] = 1.85
      elif molecule.elements[i] == 'S': parm.parm_data['RADII'][i] = 1.8
      elif molecule.elements[i] == 'Cl': parm.parm_data['RADII'][i] = 1.5
      else: parm.parm_data['RADII'][i] = 1.5
   
   parm.parm_data['RADIUS_SET'][0] = 'amber6 modified Bondi radii (amber6)'
   _screen1(parm, molecule)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def mbondi(parm):
   """ Sets the mbondi radii """
   molecule = parm.ToMolecule()
   for i in range(parm.pointers['NATOM']):
      # Radius of H atom depends on element it is bonded to
      if molecule.elements[i] == 'H':
         if molecule.elements[molecule.bonds[i][0]] in ['C', 'N']:
            parm.parm_data['RADII'][i] = 1.3
         elif molecule.elements[molecule.bonds[i][0]] in ['O', 'S']:
            parm.parm_data['RADII'][i] = 0.8
         else:
            parm.parm_data['RADII'][i] = 1.2
      # Radius of C atom depends on what type it is
      elif molecule.elements[i] == 'C':
         if parm.parm_data['AMBER_ATOM_TYPE'][i][:2] in ['C1', 'C2', 'C3']:
            parm.parm_data['RADII'][i] = 2.2
         else:
            parm.parm_data['RADII'][i] = 1.7
      # All other elements have fixed radii for all types/partners
      elif molecule.elements[i] == 'N': parm.parm_data['RADII'][i] = 1.55
      elif molecule.elements[i] == 'O': parm.parm_data['RADII'][i] = 1.5
      elif molecule.elements[i] == 'F': parm.parm_data['RADII'][i] = 1.5
      elif molecule.elements[i] == 'Si': parm.parm_data['RADII'][i] = 2.1
      elif molecule.elements[i] == 'P': parm.parm_data['RADII'][i] = 1.85
      elif molecule.elements[i] == 'S': parm.parm_data['RADII'][i] = 1.8
      elif molecule.elements[i] == 'Cl': parm.parm_data['RADII'][i] = 1.5
      else: parm.parm_data['RADII'][i] = 1.5
   
   parm.parm_data['RADIUS_SET'][0] = 'modified Bondi radii (mbondi)'
   _screen1(parm, molecule)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def mbondi2(parm):
   """ Sets the mbondi2 radii """
   molecule = parm.ToMolecule()
   for i in range(parm.pointers['NATOM']):
      # Radius of H atom depends on element it is bonded to
      if molecule.elements[i] == 'H':
         if molecule.elements[molecule.bonds[i][0]] in ['N']:
            parm.parm_data['RADII'][i] = 1.3
         else:
            parm.parm_data['RADII'][i] = 1.2
      # Radius of C atom depends on what type it is
      elif molecule.elements[i] == 'C':
         if parm.parm_data['AMBER_ATOM_TYPE'][i][:2] in ['C1', 'C2', 'C3']:
            parm.parm_data['RADII'][i] = 2.2
         else:
            parm.parm_data['RADII'][i] = 1.7
      # All other elements have fixed radii for all types/partners
      elif molecule.elements[i] == 'N': parm.parm_data['RADII'][i] = 1.55
      elif molecule.elements[i] == 'O': parm.parm_data['RADII'][i] = 1.5
      elif molecule.elements[i] == 'F': parm.parm_data['RADII'][i] = 1.5
      elif molecule.elements[i] == 'Si': parm.parm_data['RADII'][i] = 2.1
      elif molecule.elements[i] == 'P': parm.parm_data['RADII'][i] = 1.85
      elif molecule.elements[i] == 'S': parm.parm_data['RADII'][i] = 1.8
      elif molecule.elements[i] == 'Cl': parm.parm_data['RADII'][i] = 1.5
      else: parm.parm_data['RADII'][i] = 1.5

   parm.parm_data['RADIUS_SET'][0] = 'H(N)-modified Bondi radii (mbondi2)'
   _screen1(parm, molecule)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def mbondi3(parm):
   """ Sets mbondi3 radii """
   mbondi2(parm) # start from mbondi2 radii
   molecule = parm.ToMolecule()
   for i in range(parm.pointers['NATOM']):
      # Adjust OE (GLU), OD (ASP), and HH/HE (ARG)
      if parm.parm_data['RESIDUE_LABEL'][molecule.residue_container[i]] in [
                                                        'GLU','ASP','GL4','AS4']:
         if parm.parm_data['ATOM_NAME'][i].startswith('OE') or \
            parm.parm_data['ATOM_NAME'][i].startswith('OD'):
            parm.parm_data['RADII'][i] = 1.4
      elif parm.parm_data['RESIDUE_LABEL'][molecule.residue_container[i]] == 'ARG':
         if parm.parm_data['ATOM_NAME'][i].startswith('HH') or \
            parm.parm_data['ATOM_NAME'][i].startswith('HE'):
            parm.parm_data['RADII'][i] = 1.17
      # Adjust carboxylate O radii on C-Termini. Don't just do the end residue,
      # since we can have C-termini in the middle as well (i.e. 2-chain dimers)
      if parm.parm_data['ATOM_NAME'][i] == 'OXT':
         parm.parm_data['RADII'][i] = 1.4
         parm.parm_data['RADII'][i-1] = 1.4

   parm.parm_data['RADIUS_SET'][0] = \
            'ArgH and AspGluO modified Bondi2 radii (mbondi3)'
   _screen1(parm, molecule)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _screen1(parm, molecule):
   """ Applies the first set of screening parameters found in tleap source """
   for i in range(parm.pointers['NATOM']):
      if molecule.elements[i] == 'H':
         parm.parm_data['SCREEN'][i] = 0.85
      elif molecule.elements[i] == 'C':
         parm.parm_data['SCREEN'][i] = 0.72
      elif molecule.elements[i] == 'N':
         parm.parm_data['SCREEN'][i] = 0.79
      elif molecule.elements[i] == 'O':
         parm.parm_data['SCREEN'][i] = 0.85
      elif molecule.elements[i] == 'F':
         parm.parm_data['SCREEN'][i] = 0.88
      elif molecule.elements[i] == 'P':
         parm.parm_data['SCREEN'][i] = 0.86
      elif molecule.elements[i] == 'S':
         parm.parm_data['SCREEN'][i] = 0.96
      else:
         parm.parm_data['SCREEN'][i] = 0.8

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _screen2(parm, molecule):
   """
   Applies the second set of screening parameters found in tleap source
   (unused as far as I can tell)
   """
   for i in range(parm.pointers['NATOM']):
      if molecule.elements[i] == 'H':
         parm.parm_data['SCREEN'][i] = 0.8461
      elif molecule.elements[i] == 'C':
         parm.parm_data['SCREEN'][i] = 0.9615
      elif molecule.elements[i] == 'N':
         parm.parm_data['SCREEN'][i] = 0.9343
      elif molecule.elements[i] == 'O':
         parm.parm_data['SCREEN'][i] = 1.0088
      elif molecule.elements[i] == 'Na':
         parm.parm_data['SCREEN'][i] = 1.0000
      elif molecule.elements[i] == 'Mg':
         parm.parm_data['SCREEN'][i] = 1.0000
      elif molecule.elements[i] == 'P':
         parm.parm_data['SCREEN'][i] = 1.0700
      elif molecule.elements[i] == 'S':
         parm.parm_data['SCREEN'][i] = 1.1733
      else:
         parm.parm_data['SCREEN'][i] = 0.8

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _screen3(parm, molecule):
   """
   Applies the third and final set of screening parameters found in tleap source
   (unused as far as I can tell)
   """
   for i in range(parm.pointers['NATOM']):
      if molecule.elements[i] == 'H':
         parm.parm_data['SCREEN'][i] = 0.8846
      elif molecule.elements[i] == 'C':
         parm.parm_data['SCREEN'][i] = 0.9186
      elif molecule.elements[i] == 'N':
         parm.parm_data['SCREEN'][i] = 0.8733
      elif molecule.elements[i] == 'O':
         parm.parm_data['SCREEN'][i] = 0.8836
      elif molecule.elements[i] == 'Na':
         parm.parm_data['SCREEN'][i] = 1.0000
      elif molecule.elements[i] == 'Mg':
         parm.parm_data['SCREEN'][i] = 1.0000
      elif molecule.elements[i] == 'P':
         parm.parm_data['SCREEN'][i] = 0.9604
      elif molecule.elements[i] == 'S':
         parm.parm_data['SCREEN'][i] = 0.9323
      else:
         parm.parm_data['SCREEN'][i] = 0.8

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def ChRad(parm, radii_set):
   if not radii_set in ["mbondi", "mbondi2", "mbondi3", "bondi", "amber6"]:
      raise ChangeRadiiError("You must choose from mbondi, mbondi2, mbondi3, " +
                             "bondi, or amber6 radii sets!")
   
   eval("%s(parm)" % radii_set)
