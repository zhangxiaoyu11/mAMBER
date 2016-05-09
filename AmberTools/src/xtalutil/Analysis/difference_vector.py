#! /usr/bin/env python
import sys
import os
from numpy import *
from ReadAmberFiles import *
import argparse

#######################################################################################
# This script will calculate the translate vector to center trajectory on crystal
# structure.
# Arguments:
#	file1 - name of crystal rst7 file. This file must have amber prmtop atoms
#			(ie hydrogens added and same atom order) and coordinates must be
#			the same as the original crystal pdb file (centered on pdb).
#	file2 - the trajectory that will be translated. First frame of must be the
#			crystal supercell.
#	 atom -	number of atom to use to get the translation vector. I usually don't use
#			atoms from 1st or last residue. Usually I just take the first heavy atom
#			of the second residue.
# Return:
#	 Prints the move vector to file 'tmp' which can be used in bash script for ptraj.
########################################################################################

parser = argparse.ArgumentParser()
parser.add_argument("file1", help="name of file 1")
parser.add_argument("file2", help="name of file 2")
parser.add_argument("atom", help="atom number",type=int)
args = parser.parse_args()



A=rst7(args.file1)
A.coords=A.Get_Coords()
coord_pdb=A.coords[args.atom-1]

B=nc(args.file2)
B.coords=B.Get_Traj(99999)
coord_traj=B.coords[0,args.atom-1,:]


vec=coord_pdb-coord_traj
print "Translation vector is:"
print "%s\n" %vec
f=open('tmp','w')
f.write('%10.5f %10.5f %10.5f\n' %(vec[0],vec[1],vec[2]))
f.close()
