#!/bin/bash

. ../MasterTest.sh

CleanFiles rms.in rmscorr.dat 

INPUT="-i rms.in"
cat > rms.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
strip !(:2-12@CA)
rmsavgcorr out rmscorr.dat R2-12
EOF
RunCpptraj "RmsAvgCorr"
DoTest rmscorr.dat.save rmscorr.dat
CheckTest
EndTest

exit 0
