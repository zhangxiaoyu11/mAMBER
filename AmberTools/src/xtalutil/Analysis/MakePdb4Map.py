#! /usr/bin/python
import sys
import os
from numpy import *
import Scientific.IO.NetCDF
from Scientific.IO import NetCDF as Net
import time
import RunTime

#======================================================================#
#                                                                      #
# Prepare pdb frames from a trajectory for md2map.sh.                  #
# The trick is to make sure that column 13 of the pdb file has correct #
# element names. Run it once on just one frame (set frames=1) to get   #
# output of unrecognized element names. (For example, tip4p EPW's will #
# need to be removed.                                                  #
#                                                                      #
#======================================================================#


#======================================================================#
#
# SET VARIABLES
SC_TOPO='/home/pjanowsk/Case/4lzt/RunSi/4lztSi.prmtop'
SC_TRAJ='fit_wat.nc'
frames=7426
startframe=1
offsetframe=1
angles=" 88.52 108.53 111.89 P 1           1"
#
#======================================================================#


#make directory
os.system('mkdir -p PDBData')
os.chdir('PDBData')

#get cell lenghts
ofile = Net.NetCDFFile('../fit_wat.nc', 'r')
celllen=ofile.variables['cell_lengths']

#~ # Optional
#~ # prepare the experimental pdb (using supercell without added solvent)
#~ f=open('ctraj_frame', 'w')
#~ f.write('parm ../../4lztSi.prmtop\n')
#~ f.write('trajin /home/pjanowsk/Case/pepsim/RunCase/equivalency/mergtraj_equiv.nc 1 1 1\n')
#~ f.write('trajout RunCase000001.pdb pdb\n')
#~ f.close()
#~ os.system('cpptraj -i ctraj_frame')
#~ 
#~ cell=celllen[0]
#~ os.system('sed -i \'1i CRYST1   %6.3f   %6.3f   %6.3f 116.41  95.53  93.16 P 1           2  \' RunCase000001.pdb' %(cell[0],cell[1],cell[2]))
#~ 
#~ f=open('RunCase000001.pdb','r')
#~ Main= [l.strip() for l in f.readlines()]
#~ f.close()
#~ f=open('new.pdb','w')
#~ for i in Main:
	#~ if i[0:4]=='ATOM':
		#~ if i[13] not in ['H','C','O','N']:
			#~ f.write('%s           H\n' %(i))
		#~ else:
			#~ f.write('%s           %s\n' %(i,i[13]))   
	#~ else:
		#~ f.write('%s\n' %i)
#~ f.close()
#~ os.system('mv new.pdb RunCase000001.pdb')


# prepare the simulation pdbs
start_time=time.time()
for j in range(frames):
	frame=startframe+offsetframe*j
	#REPORT TIME
	frac_complete=float(frame)/frames
	print "frame: %d" %frame
	print RunTime.main(start_time=start_time, frac_complete=frac_complete)
	#CPPTRAJ CALL TO GET PDB OF FRAME
	f=open('ctraj_frame', 'w')
	f.write('parm ../../4lztSi.prmtop\n')
	f.write('trajin ../fit_wat.nc %d %d 1\n' %(frame,frame))
	f.write('trajout %04d.pdb pdb\n' %frame)
	f.close()
	os.system('cpptraj -i ctraj_frame >tmp')
	#ADD CELL PARAMETERS
	cell=celllen[frame-1]
	os.system('sed -i \'1i CRYST1   %6.3f   %6.3f   %6.3f %s  \' %04d.pdb' %(cell[0],cell[1],cell[2],angles,frame))
	#REMOVE TIP4PEW VIRTUAL ATOM
	os.system('grep -v EPW %04d.pdb >tmp' %frame)
	os.system('mv tmp %04d.pdb' %frame)
	#ADD COLUMN 78 (ELEMENT TYPE)
	f=open('%04d.pdb' %frame,'r')
	Main= [l.strip() for l in f.readlines()]
	f.close()
	f=open('new.pdb','w')
	for i in Main:
		if i[0:4]=='ATOM':
			if i[13] not in ['H','C','O','N','S']:
				if i[12] not in ['H']:
					f.write('%s           ?\n' %(i))
					print "Unrecognized atom:"
					print "%s\n" %i
				else:
					f.write('%s           %s\n' %(i,i[12])) 					
			else:
				f.write('%s           %s\n' %(i,i[13]))   
		else:
			f.write('%s\n' %i)
	f.close()
	os.system('mv new.pdb %04d.pdb' %frame)

#CLEAN UP
os.system('rm -rf ctraj_frame new.pdb tmp')	


