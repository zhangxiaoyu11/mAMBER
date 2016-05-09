## Functions for the GLYCAM Run_tests.bash script that do the main testing
##
## Please see 00_README in this directory for more documentation.
##


##
## The central test function
##
RunParmset() { # This is run for each parmset being tested
  #
  # Make sure the testing environment is good
  #
  CheckAndSetupTestDirectory "$1"
  #
  # If we can't do the tests, bail
  #
  if [ "${AbortTest}" = "y" ] ; then
    NUMTEST=$((NUMTEST+1))
    ERR=$((ERR+1))
    continue
  fi

  #
  # Make sure we start clean
  #
  CleanMe
  #
  # If that's all you want, clean more and go to the next parameter set
  #
  if  [ "$CLEAN" -eq 1 ]  ; then
    RemoveOutputDirectories
    continue
  fi

  #
  # Announce our intentions
  #
  echo "**
Testing prmtop and inpcrd built by tleap with GLYCAM_${testdir}.
**" | tee $TEST_RESULTS

  #
  # Run each of the tleap input files in LEAPIN
  #
  for i in $(ls LEAPIN) ; do
    RunTLEAP  ${LEAPIN}/$i
  done

  #
  # Run dacdif to see if the files are good
  #
  count=0
  while [ "${count}" -lt ${#SAVEDIRS[@]} ] ; do
    #
    # Check the number of output files to be used later
    #
    chk="${OUTPUTDIRS["${count}"]}"
    outfiles=0
    for i in $(ls ${chk}) ; do
      outfiles=$((outfiles+1))
    done
    #
    # This is the main testing
    #
    chk="${SAVEDIRS["${count}"]}"
    savefiles=0
    for i in $(ls ${chk}) ; do
      DoTest "${SAVEDIRS["${count}"]}"/${i} "${OUTPUTDIRS["${count}"]}"/${i}
      savefiles=$((savefiles+1))
    done
    #
    # These are some checks
    #
    # First, make sure there were files to compare with.  If not, complain.
    if [ "${savefiles}" -eq 0 ] ; then
      ErrorReport "${SAVEDIRS["${count}"]}" "${OUTPUTDIRS["${count}"]}" "No files found in ${SAVEDIRS[${count}]}"
      NUMTEST=$((NUMTEST+1))
      ERR=$((ERR+1))
    fi
    # Now see if the number of files is the same in the save directory and the output directory.
    if [ "${savefiles}" -ne ${outfiles} ] ; then
      ErrorReport "${SAVEDIRS["${count}"]}" "${OUTPUTDIRS["${count}"]}" "Unequal number of files found in comparison directories."
      NUMTEST=$((NUMTEST+1))
      ERR=$((ERR+1))
    fi
    #
    # Increment counter to check the next save/out pair.
    #
    count=$((count+1))
  done

  #
  # Clean up after -- deleting output directories (e.g., PRMTOP and INPCRD) if the tests passed
  #
  EndTest
  if [ $PASS -eq $NUMTEST ] ; then
    CleanMe
    RemoveOutputDirectories
  else
    #
    # delete files that will only confuse things if there are errors
    #
    ClearFiles "${ErrorCLEANUPFILES}"
  fi
}
##


##
## The function that runs dacdif to see if files were generated properly
##
DoTest() {
  LOCALERR=0
  NUMTEST=$((NUMTEST+1))
  #
  # See if dacdif is defined -- if not, bail.
  #
  if [ -z $DACDIF ] ; then
    AbortReport "GLYCAM" "dacdif ($DACDIF) not found."
    LOCALERR=$((LOCALERR+1))
    exit
  fi
  #
  # Check to see if the files to be diff'd are present.
  #
  if [ ! -e $1 ] ; then
    ErrorReport $1 $2 "Standard file $1 not found."
    LOCALERR=$((LOCALERR+1))
    elif [ ! -e $2 ] ; then
      ErrorReport $1 $2 "Test output file $2 not found."
      LOCALERR=$((LOCALERR+1))
  fi
  ERR=$((ERR+LOCALERR))
  #
  # If there are no errors in this test, go ahead with dacdif
  #
  if [ $LOCALERR -eq 0 ] ; then
    #
    # Run dacdif, copying output to a local temporary file
    #
    $DACDIF -r 5.e-7 $1 $2 | tee ${LOCALLOG}
    # Copy results into the main log
    cat $LOCALLOG >> $DACDIFS
    #
    # See if the temporary file reports failure & act accordingly
    #
    OK=$(grep -c "FAILURE" ${LOCALLOG})
    if [ $OK -ne 0 ] ; then
      FAIL=$((FAIL+1))
    fi
  fi
}
##


##
## This function counts/records how many tests passed, failed or had errors
##
EndTest() {
  if [ $FAIL -gt 0 ] ; then
    echo "  $FAIL out of $NUMTEST comparisons failed." >> $TEST_RESULTS
  fi
  if [ $ERR -gt 0 ]  ; then
    echo "  $ERR out of $NUMTEST tests had errors." >> $TEST_RESULTS
  fi
  PASS=$((NUMTEST-$ERR-$FAIL))
  echo "  $PASS of $NUMTEST comparisons passed." >> $TEST_RESULTS
  echo ""
}
