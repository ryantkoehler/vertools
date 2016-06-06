#!/bin/sh -u
# 6/6/16 RTK; Should probably handle failures better (report, something?)

makefiles=`ls *.makefile`
for mfile in $makefiles; do
    echo "=========================================================="
    echo $mfile
    make -f $mfile
done

# Clean up 
rm *.o
