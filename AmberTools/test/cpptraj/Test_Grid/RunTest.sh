#!/bin/bash

. ../MasterTest.sh

CleanFiles ptraj.in out.dipole out.grid out.dx out.dx.2

INPUT="ptraj.in"

# dipole
TOP="../tz2.ortho.parm7"
cat > ptraj.in <<EOF
trajin ../tz2.ortho.nc
rms first :1-13
center :1-13 mass origin 
image origin center familiar
dipole out.dipole 20 0.5 20 0.5 20 0.5 :WAT
EOF
RunCpptraj "Dipole test"
DoTest out.dipole.save out.dipole
CheckTest

# grid 
TOP="../tz2.truncoct.parm7"
cat > ptraj.in <<EOF
trajin ../tz2.truncoct.nc
rms first :1-13
center :1-13 mass origin 
image origin center familiar
grid out.grid 20 0.5 20 0.5 20 0.5 :WAT@O  
grid out.dx 20 0.5 20 0.5 20 0.5 :WAT@O opendx
EOF
RunCpptraj "Grid test"
DoTest out.grid.save out.grid
DoTest out.dx.save out.dx
CheckTest

# grid dx read
TOP="../tz2.truncoct.parm7"
cat > ptraj.in <<EOF
trajin ../tz2.truncoct.nc 1 1
grid out.dx.2 readdx out.dx.save :1588 opendx  
EOF
RunCpptraj "OpenDX Grid read test"
DoTest out.dx.save out.dx.2
CheckTest

EndTest

exit 0
