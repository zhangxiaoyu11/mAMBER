#! /usr/bin/env python
from numpy import *
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from ReadAmberFiles import *
import argparse
from Scientific.IO import NetCDF

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--SCRestart", help="Supercell restart file with correct box info on last line")
parser.add_argument("-t", "--Trajectory", help="Supercell Trajectory")
parser.add_argument("-Title", help="Title prefix for plots")
args = parser.parse_args()



#GET CRYSTAL VOLUME
rst7file=rst7(args.SCRestart)
SCBox=rst7file.Get_Box()
crystvol=Get_volume(SCBox)

#GET TRAJECTORY VOLUMES
ofile = NetCDF.NetCDFFile(args.Trajectory, 'r')
coords=ofile.variables['coordinates']
frames=coords.shape[0]
lengths=ofile.variables['cell_lengths']
angles=ofile.variables['cell_angles']
data1=zeros((frames,2))
for frame in range(frames):
	box=UCbox=hstack((lengths[frame,:],angles[frame,:]))
	data1[frame,0]=frame
	data1[frame,1]=Get_volume(box)
savetxt('volume.dat',data1, fmt='%8.2f     %12.5f')


#CONVERT TO PERCENT
frames=len(data1[:,1])
minvol=min(data1[:,1])
maxvol=max(data1[:,1])
meanvol=mean(data1[:,1])
data1[:,1]=data1[:,1]/crystvol*100
minperc=min(data1[:,1])
maxperc=max(data1[:,1])
meanperc=mean(data1[:,1])


#SAVE STATISTICS TO FILE
f=open('volume.txt','w')
f.write('%s\n' %(os.getcwd()))
f.write('crystal volume = %10f\n' %crystvol)
f.write('min volume = %10f   %6.2f\n' %(minvol,minperc))
f.write('max volume = %10f   %6.2f\n' %(maxvol,maxperc))
f.write('mean volume = %10f   %6.2f\n' %(meanvol,meanperc))
f.close()


#PLOT VOLUME
# padding for tick labels
plt.rc(('xtick.major','xtick.minor','ytick.major','ytick.minor'), pad=15)
##thicker axes frame
plt.rc('axes',linewidth=4)
plt.rc('legend', fontsize=20)
plt.rc('mathtext',default='regular')

fig=plt.figure(figsize=(16, 12))
ax = fig.add_subplot(111)



x=[]
y=[]
for i in range(frames/10):  #this is to get average over each 10 frames (.2ns)
	x.append((data1[10*i,0])/1000)
	y.append(average(data1[10*i:10*i+10,1]))
ax.plot(x,y,'b', linewidth=4)

#~ x=[]
#~ y=[]
#~ for i in range(frames/10):  #this is to get average over each 10 frames (.2ns)
	#~ x.append((data2[10*i,0])/1000)
	#~ y.append(average(data2[10*i:10*i+10,1]))
#~ ax.plot(x,y,'k', linewidth=4)

#~ ax.plot(data3[:,0]/1000,data3[:,1],'y', linewidth=4)
#~ ax.plot(data3[:,0]/1000,data3[:,2],'k', linewidth=4)
#~ ax.plot(data3[:,0]/1000,data3[:,3],'b', linewidth=4)

for label in ax.xaxis.get_ticklabels():
	label.set_fontsize(24)
for label in ax.yaxis.get_ticklabels():
	label.set_fontsize(24)
plt.title(args.Title,fontsize=28)
plt.xlabel('Time (ns)',fontsize=28, labelpad=10)
#ax.set_xticklabels([0,1.0,2.0,3.0,4.0,5.0])
plt.ylabel(r'Volume ($\%$)',fontsize=28, labelpad=10)
plt.ylim((99.5,100.5))
#plt.xlim(xmax=500)	#modify to trajectory length

ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
for line in ax.get_xticklines() + ax.get_yticklines():
	line.set_markeredgewidth(4)
	line.set_markersize(10)

#~ ax.legend(["lattice me", "lattice dave"],bbox_to_anchor=(0, 0, .95, .2))


#plt.show()
plt.savefig('volume.png') 
