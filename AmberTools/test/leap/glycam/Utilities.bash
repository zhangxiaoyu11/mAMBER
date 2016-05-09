## Utility functions for the GLYCAM Run_tests.bash script
##
## Please see 00_README in this directory for more documentation.
##


##
## A little function to run tleap
##
RunTLEAP() {
  $TLEAP -f $1 >> $OUTPUT 2>>$ERROR
}
##


##
## Functions to facilitate error reporting
##
ErrorReport() {
 echo "During the test of files ${1} and ${2}, there was an
Error:  ${3}
==============================================================" | tee -a $TEST_RESULTS
}
AbortReport() {
 echo "The ${1} test must be aborted because:
Error:  ${2}
==============================================================" | tee -a $TEST_RESULTS
}
##


##
## Functions for removing things
##
CleanMe() { # Cleans a parameter set directory
  ClearFiles "${CLEANUPFILES}"
  for i in "${OUTPUTDIRS[@]}" ; do
    if [ -d ${i} ] ; then
      ClearDirectories ${i}
    fi
  done
}
ClearFiles() { # removes all files in a list
  for i in $1 ; do
  if [ -e $i ] ; then
    rm $i
  fi
done
}
RemoveOutputDirectories() {
  for i in "${OUTPUTDIRS[@]}" ; do
    if [ -d ${i} ] ; then
      ClearDirectories ${i}
      rmdir ${i}
    fi
  done
}
ClearDirectories() { # clear all contents in a directory
  for i in ${1} ; do
    if [ -d ${i} ] ; then
      (cd ${i} && CleanDirectory );
    fi
  done
}
CleanDirectory() { # remove all files in the current directory
  for i in $(ls) ; do
    if [ -e $i ] ; then
      rm $i
    fi
  done
}
##


##
## Check the test directories and contents.  Make the output directories.
##
CheckAndSetupTestDirectory() {
  #
  # Test that the required directories are present.  If present, make sure they are directories.
  #
  for chk in "${REQUIREDDIRECTORIES[@]}" ; do
    if [ -e ${chk} ] && [ ! -d ${chk} ] ; then
      AbortReport "${testdir}" "A non-directory entity named ${chk} exists in $(pwd)"
      AbortTest="y"
    fi
    if [ ! -e ${chk} ] ; then
      AbortReport "${testdir}" "The required directory ${chk} does not exist in $(pwd)"
      AbortTest="y"
    fi
  done
  for chk in "${OUTPUTDIRS[@]}" ; do
    if [ -e ${chk} ] && [ ! -d ${chk} ] ; then
      AbortReport "${testdir}" "A non-directory entity named ${chk} exists in $(pwd)"
      AbortTest="y"
    fi
  done
  #
  # Make output directories that don't already exist.  Abort test if that isn't possible.
  #
  for chk in "${OUTPUTDIRS[@]}" ; do
    LocalAbortTest="n"
    if [ ! -d ${chk} ] && [ "${CLEAN}" -eq 0 ] ; then
      mkdir ${chk} || LocalAbortTest="y"
      if [ "${LocalAbortTest}" = "y" ] ; then
        AbortReport $1 "Cannot create ${chk} Directory."
        AbortTest="y"
        LocalAbortTest="n"
      fi
    fi
  done
  #
  # Make sure that the LEAPIN directory contains at least one file
  #
  numLEAPIN=0
  if [ -d "${LEAPIN}" ] && [ "${CLEAN}" = 0 ] ; then
    for chk in $(ls ${LEAPIN}) ; do
      numLEAPIN=$((numLEAPIN+1))
    done
    if [ "${numLEAPIN}" -eq 0 ] ; then
      AbortReport $1 "There are no leap input files to test."
      AbortTest="y"
    fi
  fi
}
