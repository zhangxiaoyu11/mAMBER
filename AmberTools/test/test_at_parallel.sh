#!/bin/sh

date_string=`date +%Y-%m-%d_%H-%M-%S`
logdir="${AMBERHOME}/logs/test_at_parallel"
logprefix="${logdir}/${date_string}"
logfile="${logprefix}.log"
difffile="${logprefix}.diff"

if [ -z "${DO_PARALLEL}" ] ; then
	echo "Error: The environment variable DO_PARALLEL is not set. Exiting."
	exit
fi

if [ -z "${AMBERHOME}" ] ; then
   dir=`dirname \`pwd\``
	echo "Error: AMBERHOME is not set! Set AMBERHOME to `dirname $dir` and re-run the tests."
	exit
fi

mkdir -p ${logdir}
(make -k -f Makefile test.parallel2 2>&1) | tee ${logfile}

passed_count=`grep PASS ${logfile} | wc -l`
questionable_count=`grep "FAILURE:" ${logfile} | wc -l`
error_count=`grep "\<[Ee]rror\>" ${logfile} | wc -l`

echo "${passed_count} file comparisons passed" | tee -a ${logfile}
echo "${questionable_count} file comparisons failed" | tee -a ${logfile}
echo "${error_count} tests experienced errors" | tee -a ${logfile}

echo "Test log file saved as ${logfile}" | tee -a ${logfile}

if [ -f TEST_FAILURES.diff ]; then
   mv TEST_FAILURES.diff ${difffile}
   echo "Test diffs file saved as ${difffile}" | tee -a ${logfile}
else
   echo "No test diffs to save!" | tee -a ${logfile}
fi

# save summary for later reporting:
tail -5 ${logfile} > ${logdir}/at_summary
