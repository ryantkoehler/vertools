#!/usr/bin/env bash
#
# 11/4/16 RTK; Simply download and build 
# 2022-07-15 RTK; Update (target version, script somewhat)

startdir=$(pwd)
temp=$(mktemp)

# Could be whatever version?
# 2022-07-15 version
url="https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_5_x/ViennaRNA-2.5.1.tar.gz"

# tarfile should be last part of url
tarfile=$(echo $url | tr "/" " " | awk '{print $NF}')


# Check if already have tar.gz 
find -maxdepth 1 -name "ViennaRNA*tar.gz" > $temp
n_targz=$(wc $temp | awk '{print $1}')
if [[ n_targz -gt 0 ]]; then
    echo "Existing package file found"
    cat $temp
    exit 1
fi
 
echo " "
echo "Getting $url"
echo " "
wget $url
if [ ! -f $tarfile ]; then
    echo "No $tarfile = no build"
    exit 1
fi

# Extract and get vienna subdir
tar -xvzf $tarfile
subdir=$(find -maxdepth 1 -type d -name "ViennaRNA*") 
if [ ! -d $subdir ]; then
    echo "No vienna subdir from $tarfile"
    exit 1
fi

cd $subdir
outdir=$(pwd)

# Put package in local dir via "--prefix"; 
# The "--disable-openmp" is needed to build without openmp 
./configure  --prefix $outdir --disable-openmp        

echo " "
echo " "
make

# Find library and copy it where we started
cd $startdir
libname="libRNA.a"
libfile=$(find -name $libname)
if [ ! -f $libfile ]; then
    echo "No vienna library found: $libname"
    exit 1
fi

cp $libfile $startdir


echo " "
echo "Vienna library file: $libfile"
echo " "

