#!/bin/sh

# DBC454 test


if [ ! `command -v mpirun` ];   then echo "ERROR: mpi is not installed (cannot find mpirun)"; exit; fi
if [ ! -f dbc454 ];	            then echo "ERROR: package incomplete: dbc454 executable is missing"; exit; fi
if [ ! -f LICENSE_gpl-3.0.txt ];then echo "ERROR: package incomplete: LICENSE file is missing"; exit; fi
if [ ! -f LICENSE_gpl-2.0.txt ];then echo "ERROR: package incomplete: LICENSE file is missing"; exit; fi
if [ ! -f README.txt ];         then echo "ERROR: package incomplete: README file is missing"; exit; fi
if [ ! -f test.fa ];            then echo "ERROR: package incomplete: test.fa file is missing"; exit; fi
if [ ! -f unit_test.expected ]; then echo "ERROR: package incomplete: unit_test.expected file is missing"; exit; fi
if [ -f unit_test.clusters ];   then rm unit_test.clusters ; fi
if [ -f unit_test.stde ];       then rm unit_test.stde ; fi
if [ -f unit_test.stdo ];       then rm unit_test.stdo ; fi

echo "Running dbc454 test (it can take a few minutes)"
mpirun -n 8 ./dbc454 -i test.fa  -n 100 -c1 -o unit_test.clusters 2> unit_test.stde

echo "Parsing results"
grep "^DBV" unit_test.stde | cut -d, -f2- > unit_test.observed

echo "Checking results"
LCNT=`wc -l unit_test.clusters`
LCNT=`echo $LCNT | cut -d\   -f1`

if [ $LCNT == 100001 ]; then
	echo "PASSED: file 'unit_test.clusters' is complete"
	rm unit_test.clusters
else
	echo "FAILED: file 'unit_test.clusters' has $LCNT lines instead of 100001 lines, as expected."
fi

ERR=`diff unit_test.observed unit_test.expected | wc -l`
if [ $ERR == 0 ]; then
	echo "PASSED: Clustering Process was Successfull"
	rm unit_test.stde
	rm unit_test.observed
else
	echo "FAILED: Validation failed; see differences with expected results:"
	diff unit_test.observed unit_test.expected
        cat unit_test.stde
fi


