## Setup information for the GLYCAM Run_tests.bash script
##
## Please see 00_README in this directory for more documentation.
##


##
## Some locations and files that are always needed
##
# These are the directories for each parameter set to be tested
TESTDIRECTORIES=(06h 06EPb)
##
##########
## !!!!!!
## The following two arrays must match in number and order
 # directories for saving new output
 OUTPUTDIRS=(PRMTOP INPCRD)
 # the standard output should be in these directories
 SAVEDIRS=(save.TOP save.CRD)
## !!!!!!
##########
# The leap input files should be in this directory
LEAPIN="LEAPIN"
# The tests can not run unless these directories are present
REQUIREDDIRECTORIES=(${SAVEDIRS[@]} ${LEAPIN})
# These files can be removed after successful tests
CLEANUPFILES="TEST_RESULTS TLEAP_ERROR TLEAP_OUT DACDIF_ERROR DACDIF_OUTPUT LocalLog leap.log"
# These files should be removed even after unsuccessful tests
ErrorCLEANUPFILES="LocalLog leap.log"
##


##
## Check if we only need to clean
##
CLEAN=0
if [ "$1" = "clean" ] ; then
  CLEAN=1
fi
##


##
## If this run is not just for cleaning, then perform some basic setup
##
if [ $CLEAN -eq  0 ] ; then  # Do some other setup for the tests
  #
  # Set names for test ouput files
  #
  OUTPUT="TLEAP_OUT"  # standard output from the tleap runs
  ERROR="TLEAP_ERROR" # standard error from the tleap runs
  TEST_RESULTS="TEST_RESULTS" # save various messages from the testing here
  DACDIFS="DACDIF_OUTPUT" # standard output from the dacdif runs
  DACDIFERRORS="DACDIF_ERROR" # standard error from the dacdif runs
  LOCALLOG="LocalLog" # a temporary file containing the current dacdif output
  #
  # Setup things that depend on $AMBERHOME
  #
  if [ -z $AMBERHOME ] ; then
    echo "Error: The GLYCAM tleap tests require AMBERHOME to be defined."
    exit 1
  fi
  DACDIF="$AMBERHOME/AmberTools/test/dacdif"
  TLEAP="$AMBERHOME/bin/tleap"
  if [ ! -x $DACDIF ] ; then
    echo "$DACDIF not found."
    exit 1
  fi
  if [ ! -x $TLEAP ] ; then
    echo "$TLEAP not found."
    exit 1
  fi
fi
