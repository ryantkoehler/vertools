4/4/14 RTK
6/2/16 RTK; update v2.2.5
11/4/16 RTK; update v2.3 and download, script commands

----------------------------------------------------------------------------
Script to download and build Vienna library.

From top vertools directory:

    scripts/get_build_Vienna.sh

Note: The above script downloads Vienna and builds it all in a subdir; 
The only file that is needed for vertools is the library: "libRNA.a"
The script does not remove the download / build Vienna package.


----------------------------------------------------------------------------
Steps to set things up manually; The above script does all of this.

Download ViennaRNA source package from:

    https://www.tbi.univie.ac.at/RNA/

... Or ...

    wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_3_x/ViennaRNA-2.3.0.tar.gz

Unpack and build via supplied configure script; Put package in local
dir (here) via "--prefix"; The "--disable-openmp" is needed to build
without openmp (DNA software?)

Untar, cd in and configure

    tar -xvzf ViennaRNA-2.3.tar.gz
    cd ViennaRNA-2.3
    ./configure  --prefix $progs/ViennaFoldLib/ --disable-openmp        
    make
    make clean

Copy to code / vertools dir

    find . -name "libRNA.a" -exec cp {} $cdir \;

    find . -name "libRNA.a" -exec cp {} $gitdir/vertools \;

