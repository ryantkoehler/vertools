#!/bin/sh -eu
# 11/4/16 RTK; Simply download and build 

startdir=`pwd`

# Could be whatever version?
tarfile="ViennaRNA-2.3.0.tar.gz"
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_3_x/$tarfile

if [ ! -f $tarfile ]; then
    echo "No $tarfile = no build"
    exit 1
fi

# Extract and get vienna subdir
tar -xvzf $tarfile
rm $tarfile
subdir=`find . -type d | grep ViennaRNA`
if [ ! -f $subdir ]; then
    echo "No vienna subdir from $tarfile"
    exit 1
fi

cd $subdir
outdir=`pwd`

# Put package in local dir via "--prefix"; 
# The "--disable-openmp" is needed to build without openmp 
./configure  --prefix $outdir --disable-openmp        
make
make clean

# Find library and copy it where we started
cd $startdir
libname="libRNA.a"
find . -name $libname -exec cp {} $startdir \;
if [ ! -f $libname ]; then
    echo "No vienna library found: $libname"
    exit 1
fi

echo " "
echo "Vienna library file: $libname"
echo " "

