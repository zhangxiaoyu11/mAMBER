# This should be sourced from Test run scripts

# If arg is Key=Value, separate into Key and Value
ParseArg() {
  KEY=`echo "$1" | awk 'BEGIN{FS = "=";}{print $1;}'`
  VALUE=`echo "$1" | awk 'BEGIN{FS = "=";}{print $2;}'`
  if [[ $VALUE = $KEY ]] ; then
    VALUE=""
  fi
}

# DoTest(): Compare File1 to File2, print an error if they differ.
#           Args 3 and 4 can be used to pass an option to diff
DoTest() {
  # AmberTools - Use dacdif if defined.
  if [[ ! -z $DACDIF ]] ; then
    $DACDIF $1 $2
  else
    if [[ $NOTEST -eq 0 ]] ; then
      ((NUMTEST++))
      if [[ ! -e $1 ]] ; then
        echo "  $1 not found." >> $TEST_RESULTS
        echo "  $1 not found." >> $TEST_ERROR
        ((ERR++))
      elif [[ ! -e $2 ]] ; then
        echo "  $2 not found." >> $TEST_RESULTS
        echo "  $2 not found." >> $TEST_ERROR
        ((ERR++))
      elif [[ `diff $1 $2 $3 $4 | wc -l` -gt 0 ]] ; then
        #echo "------------------------------------------------------------" >> $TEST_ERROR
        echo "  $1 $2 are different." >> $TEST_RESULTS
        #diff $1 $2 $3 $4 >> $TEST_RESULTS 2>&1
        echo "  $1 $2 are different." >> $TEST_ERROR
        diff $1 $2 $3 $4 >> $TEST_ERROR 2>&1
        ((ERR++))
      else
        echo "  $2 OK." >> $TEST_RESULTS
      fi
    fi
  fi
}

# CheckTest(): Report if the error counter is greater than 0.
CheckTest() {
  # AmberTools - Dont do this if dacdif defined.
  if [[ -z $DACDIF ]] ; then
    if [[ $ERR -gt 0 ]] ; then
      echo "  $ERR comparisons failed so far."
      #echo "  $ERR out of $NUMTEST comparisons failed."
      #echo "  $ERR out of $NUMTEST comparisons failed." >> $TEST_RESULTS
      #echo "---------------------------------------------------------"
      #echo "---------------------------------------------------------" >> $TEST_RESULTS
      #exit 1
    #else
    # echo "  $NUMTEST comparisons ok."
    # echo "  $NUMTEST comparisons ok." >> $TEST_RESULTS
    fi
  fi
}

# RunCpptraj(): Run cpptraj with the given options. Start and stop MPI if requested.
RunCpptraj() {
  # If only cleaning requested no run needed, exit now
  if [[ $CLEAN -eq 1 ]] ; then
    exit 0
  fi
  echo ""
  echo "  CPPTRAJ: $1"
  if [[ -z $DACDIF ]] ; then
    echo "  CPPTRAJ: $1" >> $TEST_RESULTS
  fi
  #$STARTMPI
  if [[ ! -z $DEBUG ]] ; then
    echo "$TIME $DO_PARALLEL $VALGRIND $CPPTRAJ $DEBUG $TOP $INPUT >> $OUTPUT 2>>$ERROR"
  fi
  $TIME $DO_PARALLEL $VALGRIND $CPPTRAJ $DEBUG $TOP $INPUT >> $OUTPUT 2>>$ERROR
  #$STOPMPI
}

# EndTest(): Called at the end of every test script if no errors found.
EndTest() {
  # AmberTools - Dont do this if dacdif defined.
  if [[ -z $DACDIF ]] ; then
    if [[ $ERR -gt 0 ]] ; then
      echo "  $ERR out of $NUMTEST comparisons failed."
      echo "  $ERR out of $NUMTEST comparisons failed." >> $TEST_RESULTS
      echo "  $ERR out of $NUMTEST comparisons failed." >> $TEST_ERROR
    else 
      echo "All $NUMTEST comparisons passed." 
      echo "All $NUMTEST comparisons passed." >> $TEST_RESULTS 
    fi
    echo ""
    if [[ ! -z $VALGRIND ]] ; then
      echo "Valgrind summary:"
      grep ERROR $ERROR
      grep heap $ERROR
      grep LEAK $ERROR
      echo ""
      echo "Valgrind summary:" >> $TEST_RESULTS
      grep ERROR $ERROR >> $TEST_RESULTS
      grep heap $ERROR >> $TEST_RESULTS
      grep LEAK $ERROR >> $TEST_RESULTS

    fi
    if [[ $PROFILE -eq 1 ]] ; then
      if [[ -e gmon.out ]] ; then
        gprof $CPPTRAJ > profiledata.txt
      fi
    fi
    #echo "---------------------------------------------------------"
    #echo "---------------------------------------------------------" >> $TEST_RESULTS
  fi
}

# CleanFiles(): For every arg passed to the function, check for the file and rm it
CleanFiles() {
  while [[ ! -z $1 ]] ; do
    #for RMFILE in `find . -name "$1"` ; do
    if [[ -e $1 ]] ; then
      #echo "  Cleaning $1"
      rm $1
    fi
    #done
    shift
  done
  # If only cleaning requested no run needed, exit now
  if [[ $CLEAN -eq 1 ]] ; then
    exit 0
  fi
}

# Library Checks - Tests that depend on certain libraries like Zlib can run
# these to make sure cpptraj was compiled with that library - exit gracefully
# if not.
# Should not be called if CLEAN==1, CleanFiles should always be called first.
CheckZlib() {
  if [[ -z $ZLIB ]] ; then
    echo "This test requires zlib. Cpptraj was compiled without zlib support."
    echo "Skipping test."
    #echo "---------------------------------------------------------"
    exit 0
  fi
}

CheckBzlib() {
  if [[ -z $BZLIB ]] ; then
    echo "This test requires bzlib. Cpptraj was compiled without bzlib support."
    echo "Skipping test."
    #echo "---------------------------------------------------------"
    exit 0
  fi
}

CheckNetcdf() {
  if [[ -z $NETCDFLIB ]] ; then
    echo "This test requires Netcdf. Cpptraj was compiled without Netcdf support."
    echo "Skipping test."
    #echo "---------------------------------------------------------"
    exit 0
  fi
}

CheckPtrajAnalyze() {
  if [[ ! -z $NOMATHLIB ]] ; then
    echo "This test requires LAPACK/ARPACK/BLAS routines from AmberTools."
    echo "Cpptraj was compiled with -DNO_MATHLIB. Skipping test."
    #echo "---------------------------------------------------------"
    exit 0
  fi
}
#==============================================================================
# If the first argument is "clean" then no set-up is required. Script will
# exit when either CleanFiles or RunCpptraj is called from sourcing script.
CLEAN=0
if [[ $1 = "clean" ]] ; then
  CLEAN=1
fi

# If not cleaning, determine which binary to test. By default we are
# running tests with AmberTools, in which case AMBERHOME/bin/cpptraj 
# will be tested. If AMBERHOME is not defined or if standalone is
# specified then CPPTRAJHOME/bin/cpptraj will be tested. 
DACDIF=""
if [[ $CLEAN -eq 0 ]] ; then
  # Option defaults
  NODACDIF=0
  STANDALONE=0
  CPPTRAJ=""
  TIME=""
  VALGRIND=""
  STARTMPI=""
  STOPMPI=""
  NP=1
  MPI=0
  TOP=
  INPUT=
  NOTEST=0
  OUTPUT="test.out"
  PROFILE=0
  if [[ -e $OUTPUT ]] ; then
    rm $OUTPUT
  fi
  #OUTPUT="/dev/stdout"
  ERROR="/dev/stderr"
  DEBUG=""
  while [[ ! -z $1 ]] ; do
    case "$1" in
      "stdout" ) OUTPUT="/dev/stdout" ;;
      "nodacdif" ) NODACDIF=1 ;;
      "standalone" ) STANDALONE=1 ;;
      "vg"  ) 
        VG=`which valgrind`
        if [[ -z $VG ]] ; then
          echo "vg: Valgrind not found."
          echo "    Make sure valgrind is installed and in your PATH"
          exit 1
        fi
        echo "  Using Valgrind."
        VALGRIND="valgrind --tool=memcheck --leak-check=yes --show-reachable=yes" 
        ERROR="valgrind.out"
        if [[ -e $ERROR ]] ; then
          rm $ERROR
        fi
      ;;
      "vgh" )
        VG=`which valgrind`
        if [[ -z $VG ]] ; then
          echo "vg: Valgrind not found."
          echo "    Make sure valgrind is installed and in your PATH"
          exit 1
        fi
        echo "  Using Valgrind."
        VALGRIND="valgrind --tool=helgrind"
        ERROR="valgrind.out"
        if [[ -e $ERROR ]] ; then
          rm $ERROR
        fi
        ;;
      "mpi" ) MPI=1 ;;
      "time") TIME="time" ;;
      "np"  ) 
        shift
        NP=$1
      ;;
      "-i" )
        shift
        INPUT=$1
        echo "Using input file: $INPUT"
        NOTEST=1
      ;;
      "-p" )
        shift
        TOP="$TOP -p $1"
        echo "Using top file: $1"
        NOTEST=1
      ;;
      "notest" )
        echo "End of run tests will be skipped."
        NOTEST=1
      ;;
      "-d"    ) DEBUG="-debug 4" ;;
      "-debug" )
        shift
        DEBUG="-debug $1"
        ;;
      "-cpptraj" )
        shift
        CPPTRAJ=$1
        echo "Using cpptraj: $CPPTRAJ"
        ;;
      "-profile" )
        PROFILE=1
        echo "Performing gnu profiling during EndTest."
        ;;
      * ) echo "Unknown opt: $1" ;;
    esac
    shift
  done

  # Check for binary. If not defined, first check AMBERHOME/bin, then
  # CPPTRAJHOME/bin, then $AMBERHOME/AmberTools/src/cpptraj/bin.
  # If using AMBERHOME/bin binary use DACDIF for test comparisons,
  # otherwise diff will be used.
  DACDIF=""
  if [[ -z $CPPTRAJ ]] ; then
    if [[ $STANDALONE -eq 0 && ! -z $AMBERHOME ]] ; then
      DACDIF=$AMBERHOME/AmberTools/test/dacdif
      if [[ ! -e $DACDIF ]] ; then
        echo "$DACDIF not found."
        exit 1
      fi
      CPPTRAJ=$AMBERHOME/bin/cpptraj
    elif [[ ! -z $CPPTRAJHOME ]] ; then
      CPPTRAJ=$CPPTRAJHOME/bin/cpptraj
    elif [[ $STANDALONE -eq 1 && ! -z $AMBERHOME ]] ; then
      CPPTRAJ=$AMBERHOME/AmberTools/src/cpptraj/bin/cpptraj
    else
      echo "Tests require CPPTRAJHOME or AMBERHOME to be defined, or specify"
      echo "cpptraj binary location with '-cpptraj <filename>'"
      exit 0
    fi
  fi
  if [[ $NODACDIF -eq 1 ]] ; then
    DACDIF=""
  fi

  # If not doing AmberTools tests, record results of each test to a file
  TEST_RESULTS=""
  if [[ -z $DACDIF ]] ; then
    TEST_RESULTS=Test_Results.dat
    TEST_ERROR=Test_Error.dat
    if [[ -e $TEST_RESULTS ]] ; then
      rm $TEST_RESULTS
    fi
    if [[ -e $TEST_ERROR ]] ; then
      rm $TEST_ERROR
    fi
  fi

  # Check for existance of binary file
  if [[ ! -e $CPPTRAJ ]] ; then
    echo "CPPTRAJ not found ($CPPTRAJ)."
    exit 1
  fi
  if [[ ! -z $DEBUG || $STANDALONE -eq 1 ]] ; then
    ls -l -t $CPPTRAJ
  fi

  # Check how cpptraj was configured to determine whether it includes 
  # netcdf, zlib, etc
  ZLIB=""
  BZLIB=""
  NETCDFLIB=""
  MPILIB=""
  NOMATHLIB=""
  DEFINES=`$CPPTRAJ --defines | grep Compiled`
  ZLIB=`echo $DEFINES | grep DHASGZ`
  BZLIB=`echo $DEFINES | grep DHASBZ2`
  NETCDFLIB=`echo $DEFINES | grep DBINTRAJ`
  MPILIB=`echo $DEFINES | grep DMPI`
  NOMATHLIB=`echo $DEFINES | grep DNO_MATHLIB`

  # Start test results file
  echo "**************************************************************"
  echo "TEST: `pwd`"
  if [[ -z $DACDIF ]] ; then 
    echo "**************************************************************" > $TEST_RESULTS
    echo "TEST: `pwd`" >> $TEST_RESULTS
  fi

  # Set up MPI environment if specified or if >1 processor requested.
  if [[ $MPI -eq 1 || $NP -gt 1 ]] ; then
    echo "  Using MPI with $NP processors."
    # Check that cpptraj was compiled mpi, but dont abort
    if [[ -z $MPILIB ]] ; then
      echo "Warning: $CPPTRAJ was not compiled with -DMPI"
    fi
    # Check that mpi is active and in path
    mpdringtest > /dev/null
    if [[ $? -ne 0 ]] ; then
      echo "  Error: No MPI daemon is running. Aborting parallel test."
      exit 1
    fi
    #STARTMPI="mpdboot -n 1"
    #STOPMPI="mpdallexit"
    DO_PARALLEL="mpiexec -n $NP"
  fi
fi # END if CLEAN==0

NUMTEST=0
ERR=0

