#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles center.in centered.crd origin.centered.crd origin.mass.centered.crd

# Test 1
CheckNetcdf
cat > center.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 1
center :1-13
trajout centered.crd
EOF
INPUT="-i center.in"
RunCpptraj "Center command test."
DoTest centered.crd.save centered.crd

cat > center.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 1
center :1-13 origin
trajout origin.centered.crd
EOF
INPUT="-i center.in"
RunCpptraj "Center origin command test."
DoTest origin.centered.crd.save origin.centered.crd

cat > center.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 1
center :1-13 origin mass
trajout origin.mass.centered.crd
EOF
INPUT="-i center.in"
RunCpptraj "Center origin mass command test."
DoTest origin.mass.centered.crd.save origin.mass.centered.crd

CheckTest

EndTest

exit 0
