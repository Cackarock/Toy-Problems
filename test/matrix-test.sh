#!/bin/bash

if [ "$1" = "" ];then
    echo "usage: $0 <output file>"
    echo "   output file - the file to save output in"
    echo "   if output file exists, this script will append to it"
    exit 0;
fi

total=0
maxpts=0
dest="$1"
TESTFILES=./input
OUTPUTFILES=./test/test-output

addPoint() {
    let "total=total+1"
    let "maxpts=maxpts+1"
}

removePoint() {
    if [ "$total" == 0 ]; then
        let "total=0"
        let "maxpts=maxpts+1"
    else 
        let "total=total-1"
        let "maxpts=maxpts+1"
    fi
}

displayFinal() {
    echo
    echo "-----$total out of $maxpts tests PASSED-----"
    if [ "$total" != "$maxpts" ]; then
        echo "Errors can be found in matrix-test-results"
    else 
        echo "Matrix test log can be found in matrix-test-results"
    fi
    echo "$total out of $maxpts tests PASSED" >> $dest
}

if [ -f $dest ];
then
    removePoint
	echo >> $dest
	echo "unit_tests.sh: output file $dest exists, appending to it" >> $dest
	echo >> $dest
else
    addPoint
fi

echo
echo "***matrix-test.sh -- running tests***"
echo

make mm 
if [ ! $? -eq 0 ];then
    removePoint
    echo "unit_tests: FAIL - make returned non-zero"  >> $dest
    exit 1
else 
    addPoint
fi



if [ ! -x "./bin/mm" ];then
    removePoint
    echo "MAKE: FAIL - no exe named mm in ./bin/" >> $dest
    exit 1
else
    addPoint
fi

echo >> $dest
echo "START: Testing Matrix-Matrix Multiply" >> $dest

# Test for invalid matrices in Cannon's multiply
./bin/mm ./input/matrix2x2.m ./input/matrix3x3.m ./test/test-output/outputMatrix.m
if [ "$?" == 0 ]; then
    removePoint
    echo "ERROR: mm accepts invalid matrices" >> $dest
else
    addPoint
    echo "---SUCCESS: mm rejects invalid matrices" >> $dest
fi

# Test for invalid args in Cannon's multiply
./bin/mm ./input/matrix1.m hello ./test/test-output/outputMatrix.m
if [ "$?" == 0 ]; then
    removePoint 
    echo "ERROR: mm accepts invalid arguments" >> $dest
else
    addPoint
    echo "---SUCCESS: mm rejects invalid arguments" >> $dest
fi

# Test for correct multiplication
./bin/mm ./input/matrix2x2.m ./input/matrix2x2.m ./test/test-output/output2x2.m
diff -i -w -B ./test/test-output/output2x2.m $TESTFILES/expected2x2.m >>diff.out
if [ "$?" == 0 ]; then 
    addPoint
    rm diff.out
    echo "---SUCCESS: mm correctly multiplies matrix x matrix" >> $dest
else
    removePoint
    echo "ERROR: mm is not calculating correctly" >> $dest
    echo "Note: The diff result is in ./diff.out" >> $dest
    exit 1
fi

echo >> $dest
#echo "CLEANING: ---" >> $dest
#make clean
#if [ "$?" == 0 ]; then 
#    addPoint
#    echo "---SUCCESS: make clean" >> $dest
#else
#    removePoint
#    echo "ERROR: make clean did not work" >> $dest
#fi

displayFinal
