#!/bin/sh

# This script is passed the BINDIR as the first argument, the PYTHON to use as
# the second, and "par" as the third if we're building in parallel

if [ $# -ne 3 -a $# -ne 2 ]; then
   # If we're not installing python
   exit 0
fi

# Now check that we didn't get "par" as our second argument -- this means that
# we're asking to build in parallel, but don't have a python
if [ "$2" = "par" ]; then
   exit 0
fi

# Create ante-MMPBSA.py
sed -e "s@PYTHONEXE@$2@g" < ante-MMPBSA.py > $1/ante-MMPBSA.py
/bin/chmod +x $1/ante-MMPBSA.py

if [ $# -eq 2 ]; then
   lastarg=$2
else
   lastarg=$3
fi

if [ "$lastarg" = "par" ]; then

   # First check and see if we already have mpi4py installed in our Python
   # and see if it works. If it does, then there's no reason to build mpi4py
   # all over again.
cat > test_mpi4py.py << EOF
#! $2
import os, sys
if os.getenv('AMBERHOME') is not None:
   sys.path.append(os.path.join(os.getenv('AMBERHOME'), 'bin'))
from mpi4py import MPI

size = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()

print "My size is %d, my rank is %d" % (size, rank)
EOF
   
   /bin/chmod +x test_mpi4py.py
   ./test_mpi4py.py > /dev/null 2>&1
   if [ $? -ne 0 ]; then
      /bin/rm -fr mpi4py-1.2.2/
      tar zxf mpi4py-1.2.2.tar.gz
      cd mpi4py-1.2.2/
      echo " Building mpi4py (this may take a while)..."
      $2 setup.py build > ../../mpi4py_install.log 2>&1
      if [ $? -gt 0 ]; then
         echo " Error in mpi4py install. Check mpi4py_install.log"
         exit 1
      else
   	   /bin/mv build/lib.*/mpi4py $1
      fi
      echo " Done building mpi4py."
      cd ..
   fi
   /bin/rm -f test_mpi4py.py

   # Create the MMPBSA.py.MPI program
   sed -e "s@PYTHONEXE@$2@g" MMPBSA.py > $1/MMPBSA.py.MPI
   /bin/chmod +x $1/MMPBSA.py.MPI

else # serial version

   # Create the MMPBSA.py serial program
   sed -e "s@PYTHONEXE@$2@g" MMPBSA.py > $1/MMPBSA.py
   /bin/chmod +x $1/MMPBSA.py
fi

/bin/cp -LR MMPBSA_mods $1
