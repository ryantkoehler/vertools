#!/bin/sh -eu
# 11/4/16 RTK; Simply find and move exectuable files

if [ $# -lt 1 ]; then
    echo "Use: <outdir>"
    echo "Finds executable files in current dir and moves them to <outdir>"
    exit 1
fi

odir=$1
if [ ! -d $odir ] || [ ! -w $odir ]; then
    echo "Probelem with output dir: $odir"
    echo "Does not exist or is not writeable"
    exit 1
fi

find . -maxdepth 1 -type f -perm -100 -exec mv {} $odir \;

