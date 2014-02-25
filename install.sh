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
cd ..
git clone https://github.com/acado/acado.git -b stable ACADOtoolkit
cd ACADOtoolkit
git checkout v1.2.1beta
mkdir build
cd build
cmake ..
make
cd ../../../
make
cd examples
make
