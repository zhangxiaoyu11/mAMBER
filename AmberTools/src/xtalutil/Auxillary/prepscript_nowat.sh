#!/bin/bash -f 
set -x
#======================================================================#

ROOT=nowat




#======================================================================#
cp ${ROOT}.pdb ${ROOT}_Orig.pdb
reduce -trim ${ROOT}.pdb >${ROOT}_NoH.pdb 2>trim.err

# Prepare coordinates for xtal simulation  AMBER
cat > rnatleap << EOF
source leaprc.ff10
loadamberparams frcmod.ionsjc_tip4pew
set default nocenter on
loadAmberPrep MPD.prepin
x = loadpdb "${ROOT}_NoH.pdb"
bond x.40.SG x.95.SG
bond x.26.SG x.84.SG
bond x.65.SG x.72.SG
bond x.58.SG x.110.SG
bond x.167.SG x.222.SG
bond x.153.SG x.211.SG
bond x.192.SG x.199.SG
bond x.185.SG x.237.SG
set x box 60
saveAmberParm x ${ROOT}.prmtop ${ROOT}.rst7
quit
EOF

tleap -f rnatleap > tleap.out
ambpdb -p ${ROOT}.prmtop < ${ROOT}.rst7 > ${ROOT}.pdb

# Set the unit cell dimensions
ChBox \
  -c ${ROOT}.rst7 \
  -o ${ROOT}.rst7 \
  -X 30.000 \
  -Y 38.270 \
  -Z 53.170 \
  -al 90.00  \
  -bt 106.00 \
  -gm 90.00
