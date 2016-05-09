#!/bin/sh

# This script is invoked via [setup.sh $BINDIR $PYTHON] which passes us
# where the programs should be installed and which python to use. If either
# is missing, we just stop

if [ $# -lt 2 ]; then
   exit 0
fi

# This means we don't have a Python
if [ $# -eq 2 -a "$2" = "par" ]; then
   exit 0
fi

# This means we're building in parallel
if [ $# -eq 3 ]; then
   par="yes"
else
   par="no"
fi

# chemistry is needed by MMPBSA.py.MPI, so we need to install the chemistry/
# package in parallel
$2 -c "from chemistry.amber import readparm" || exit 1
/bin/cp -LR chemistry $1/


if [ "$par" = "no" ]; then
   $2 -c "from ParmedTools.ParmedActions import *" || exit 1
   $2 -c "from ParmedTools import *" || exit 1

   /bin/cp -LR ParmedTools $1/
   sed -e "s@PYTHONEXE@$2@g" < parmed.py > $1/parmed.py
   sed -e "s@PYTHONEXE@$2@g" < xparmed.py > $1/xparmed.py

   /bin/chmod +x $1/parmed.py $1/xparmed.py
fi
