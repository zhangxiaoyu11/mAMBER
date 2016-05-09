#! /bin/bash
set -x

#----------------------------------------------------------------------#
#                                                                      #
# A preparation script to merge and strip of waters trajectory files   #
# before passing to FullAnalysis.sh from xtalutil/Analysis package     #
#                                                                      #
#----------------------------------------------------------------------#

#SET VARIABLES
WORKING_DIR=/home/pjanowsk/Case/4lzt/RunSi/ptraj
TRAJ_ROOT=/home/case/xtal/4lzt_xtal/mdi
TRAJ_EXT=nc
TRAJ_START=1
TRAJ_END=75
OFFSET=10
STRIP_MASK=":1669-9999"
SC_TOPO=/home/pjanowsk/Case/4lzt/RunSi/4lztSi.prmtop


########################################################################

cd ${WORKING_DIR}
rm -rf ctraj.merge.in merge_nowat.nc

# strip waters and mergetrajectory
echo -e '\n############################\nmerging trajectory\n####################'
echo "parm ${SC_TOPO}" >> ctraj.merge.in
for i in `seq ${TRAJ_START} ${TRAJ_END}`; do
	echo "trajin ${TRAJ_ROOT}${i}.${TRAJ_EXT} 1 -1 ${OFFSET}" >> ctraj.merge.in
done
echo "strip ${STRIP_MASK}" >> ctraj.merge.in
echo "center mass origin" >> ctraj.merge.in
echo "trajout merge_nowat.nc netcdf" >> ctraj.merge.in
cpptraj <ctraj.merge.in
