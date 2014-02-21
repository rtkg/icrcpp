git submodule init 
git submodule update 
cd tools/obj-1.0/
./configure --prefix=`pwd`/build
make 
make install
cd ../Qhull-2011.2/
cmake -DCMAKE_INSTALL_PREFIX=`pwd`/build
make 
make install
make clean
cd ../../
make
cd examples
make
