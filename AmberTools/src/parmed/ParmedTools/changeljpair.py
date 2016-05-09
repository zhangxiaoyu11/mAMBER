""" Changes the LJ coefficient terms only for a specific type-pair """

def ChLJPair(parm, atom_1, atom_2, rmin, eps, one_4=False):
   
   if one_4: key = 'LENNARD_JONES_14'
   else:     key = 'LENNARD_JONES'

   # Make sure that atom type 1 comes first
   if atom_1 <= atom_2:
      atm_1_idx = atom_1
      atm_2_idx = atom_2
   else:
      atm_1_idx = atom_2
      atm_2_idx = atom_1
   
   # Find the atom1 - atom2 interaction (adjusting for indexing from 0)
   term_idx = parm.parm_data['NONBONDED_PARM_INDEX'][
                   parm.pointers['NTYPES']*(atm_1_idx - 1) + atm_2_idx - 1] - 1
   
   # Now change the ACOEF and BCOEF arrays, assuming the proper amber combining
   # rules
   parm.parm_data['%s_ACOEF' % key][term_idx] = eps * rmin ** 12
   
   parm.parm_data['%s_BCOEF' % key][term_idx] = 2 * eps * rmin ** 6
