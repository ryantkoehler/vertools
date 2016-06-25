4/4/14 RTK
6/2/16 RTK; update v2.2.5

Download ViennaRNA source package from:

    https://www.tbi.univie.ac.at/RNA/


Unpack and build via supplied configure script; Put package in local
dir (here) via "--prefix"; The "--disable-openmp" is needed to build
without openmp (DNA software?)

Untar, cd in and configure

    tar -xvzf ViennaRNA-2.2.5.tar.gz
    cd ViennaRNA-2.2.5
    ./configure  --prefix $progs/ViennaFoldLib/ --disable-openmp        
    make
    make clean

Copy to code / vertools dir

    find . -name "libRNA.a" -exec cp {} $cdir \;

    find . -name "libRNA.a" -exec cp {} $gitdir/vertools \;

