#!/bin/bash 
git submodule init 
git submodule update 
cd tools/obj-1.0/
./configure --prefix=`pwd`/build
make 
make install
cd ../Qhull-2011.2/
make
make clean
cd ../../
make
cd examples
make
