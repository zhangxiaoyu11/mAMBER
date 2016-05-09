#! /usr/bin/env python
import sys
import os
from numpy import *
import commands
import argparse
from Scientific.IO import NetCDF

########################################################################

# This script uses ptraj to calculate a distance matrix first for the pdb 
# structure and than average distances for trajectory blocks. It than subtracts
# one distance matrix from the pdb distance matrix to get a difference distance
# matrix.

# Plots are made of the entire supercell distance matrix and of a single 
# asymmetric unit difference distance matrix averaged over all the copies of that
# asymmetric unit in the supercell differenc distance matrix.

# This script needs to get carefully modified for a new system. Atom selections 
# and trajectory block (start and end frames) in the ptraj scripts need to 
# be changed. By default it separates the trajectory into three blocks and
# calculates C-alpha distances.
########################################################################

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--Trajectory", help="Supercell Trajectory")
parser.add_argument("-p", "--SCTopo", help="Supercell Topology")
parser.add_argument("-u", "--UnitCells", help="number of unit cells in crystal")
parser.add_argument("-a", "--ASUs", help="number of asymmetric units per unit cell")
parser.add_argument("-tm", "--Timestep", help="Timestep for rmsd output")
args = parser.parse_args()


###########################
# SETUP                   #
###########################
print args.Trajectory, args. SCTopo, args.UnitCells, args.ASUs, args.Timestep
TRAJ=args.Trajectory
topo=args.SCTopo
uc=int(args.UnitCells)
asym=int(args.ASUs)
tstep=float(args.Timestep)


########################################################################
#                                                                      #
# CALCULATE DISTANCE MATRICES                                          #
#                                                                      #
########################################################################

# change into directory
os.system('mkdir -p distance_matrix')
os.chdir('distance_matrix')
ofile = NetCDF.NetCDFFile(TRAJ, 'r')
coords=ofile.variables['coordinates']
frames=coords.shape[0]
framesperblock=frames/3.0

# pdb distance matrix
f=open('ctraj_distmatrix', 'w')
f.write('parm %s\n' %topo)
f.write('trajin %s 1 1\n' %TRAJ)
f.write('matrix dist @CA out dist_pdb.dat\n')
f.close()
os.system('cpptraj <ctraj_distmatrix')

#pdb avg distance matrices for 20ns blocks
i=1
start=frames-framesperblock*i
end=frames-framesperblock*(i-1)
f=open('ctraj_distmatrix', 'w')
f.write('parm %s\n' %topo)
f.write('trajin %s %d %d\n' %(TRAJ,start,end))
f.write('matrix dist @CA out dist_lastthird.dat\n')
f.close()
os.system('cpptraj <ctraj_distmatrix')

i=2
start=frames-framesperblock*i
end=frames-framesperblock*(i-1)
f=open('ctraj_distmatrix', 'w')
f.write('parm %s\n' %topo)
f.write('trajin %s %d %d\n' %(TRAJ,start,end))
f.write('matrix dist @CA out dist_middlethird.dat\n')
f.close()
os.system('cpptraj <ctraj_distmatrix')

i=3
start=frames-framesperblock*i
end=frames-framesperblock*(i-1)
f=open('ctraj_distmatrix', 'w')
f.write('parm %s\n' %topo)
f.write('trajin %s %d %d\n' %(TRAJ,start,end))
f.write('matrix dist @CA out dist_firstthird.dat\n')
f.close()
os.system('cpptraj <ctraj_distmatrix')




########################################################################
#                                                                      #
# CALCULATE DIFFERENCE DISTANCE MATRICES                               #
#                                                                      #
########################################################################

# diff matrices: difference between the avg distance and the pdb distance
names=commands.getoutput('ls dist*dat').split()
for i in range(len(names)):
	if 'pdb' in names[i]:
		names.pop(i)
distpdb=genfromtxt('dist_pdb.dat')
for i in names:
	disttraj=genfromtxt(i)
	savetxt('diff'+i[i.index('_'):],(disttraj-distpdb),fmt='%8.3f')

# average per monomer difference matrix (average diff matrix over the 12 monomer)
names=commands.getoutput('ls diff*dat').split()
for i in names:
	big=genfromtxt(i)
	natoms=big.shape[0]
	residues=natoms/(uc*asym)
	assert big.shape[0] == big.shape[1]
	small=zeros((residues,residues))
	for j in range(residues):
		for k in range(residues):
			total=[]
			for unit in range(uc*asym):
					#print j+unit*residues,k+unit*residues
					total.append(big[j+unit*residues,k+unit*residues])
			small[j,k]=mean(total)
	savetxt('MonomAveg'+i[i.index('_'):],small,fmt='%8.3f')



########################################################################
#                                                                      #
# PLOT MATRICES                                                        #
#                                                                      #
########################################################################

import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import copy
import math

# This function rescales the colorbar to be centered at zero.
def cmap_center_at_zero(cmap, array):
	array_range=array.min(), array.max()
	center=0.
	if not ((array_range[0] < center) and (center < array_range[1])):
		return cmap
	center_ratio=abs(center - array_range[0]) / abs(array_range[1] - array_range[0])
	if not (0. < center_ratio) & (center_ratio < 1.):
		return cmap
	a = math.log(center_ratio) / math.log(0.5)
	if a < 0.:
		return cmap
	cdict = copy.copy(cmap._segmentdata)
	fn = lambda x : (x[0]**a, x[1], x[2])
	for key in ('red','green','blue'):
		cdict[key] = map(fn, cdict[key])
		cdict[key].sort()
		assert (cdict[key][0]<0 or cdict[key][-1]>1), "Resulting indices extend out of the [0, 1] segment."
	return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)

# Supercell difference distance matrix
names=commands.getoutput('ls MonomAveg*dat').split()
for x in names:
	majorLocator   = MultipleLocator(20)
	fig=plt.figure()
	plt.suptitle(x)
	dist=genfromtxt(x)
	ax = fig.add_subplot(111)

	im=plt.imshow(dist,cmap=cmap_center_at_zero(plt.cm.seismic, dist), origin='lower', interpolation='nearest')
	#ax.yaxis.set_major_locator(majorLocator)
	plt.colorbar()
	#~ plt.show()
	plt.savefig('../'+x+'.png')

# Single asymmetric unit difference distance matrix
names=commands.getoutput('ls diff*dat').split()
for x in names:
	majorLocator   = MultipleLocator(61)
	fig=plt.figure()
	plt.suptitle(x)
	dist=genfromtxt(x)
	ax = fig.add_subplot(111)

	im=plt.imshow(dist,cmap=cmap_center_at_zero(plt.cm.seismic, dist), origin='lower', interpolation='nearest')
	plt.colorbar()
	plt.grid(True,ls="-")
	ax.set_xticks(arange(30,730,61),minor=True)
	ax.set_xticks(arange(0,730,61))
	ax.set_xticklabels(arange(1,13),minor=True)
	ax.set_xticklabels([])
	ax.set_yticks(arange(30,730,61),minor=True)
	ax.set_yticks(arange(0,730,61))
	ax.set_yticklabels(arange(1,13),minor=True)
	ax.set_yticklabels([])
	for tick in ax.xaxis.get_minor_ticks():
		tick.tick1line.set_markersize(0)
		tick.tick2line.set_markersize(0)
		tick.label1.set_horizontalalignment('center')
	for tick in ax.yaxis.get_minor_ticks():
		tick.tick1line.set_markersize(0)
		tick.tick2line.set_markersize(0)
		tick.label1.set_verticalalignment('center')
	#~ plt.show()
	plt.savefig('../'+x+'.png')
