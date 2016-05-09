#!/bin/bash
# This is the main script for the GLYCAM tests
# It was copied from the cpptraj directory, but modified heavily,
# Starting with BLFoley on 20120229
###
### Please see 00_README in this directory for more documentation.
###
. ./Setup.bash
. ./TestFunctions.bash
. ./Utilities.bash
for testdir in "${TESTDIRECTORIES[@]}" ; do
  NUMTEST=0
  ERR=0
  FAIL=0
  PASS=0
  AbortTest="n"
 ( cd ${testdir} && RunParmset ${testdir} )
done

exit 0
