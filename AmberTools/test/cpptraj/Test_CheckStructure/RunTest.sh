#!/bin/bash

. ../MasterTest.sh

CleanFiles check.in report.dat

INPUT="-i check.in"
cat > check.in <<EOF
parm ../tz2.parm7
trajin tz2.stretched.pdb
check reportfile report.dat offset 0.7
EOF
RunCpptraj "Structure Check"
DoTest report.dat.save report.dat
CheckTest
EndTest

exit 0
