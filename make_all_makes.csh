#!/bin/csh -f
# 11/26/14 RTK

set makefiles = `ls *.makefile`
foreach mfile ( $makefiles ) 
    echo "=========================================================="
    echo $mfile
    make -f $mfile
end
