IFLAGS  = [ '-I. -I../tools/eigen-eigen-3.0.2/ -I../tools/Qhull-2011.2/src/libqhullcpp ' ];

CPPFLAGS  = [IFLAGS, '-D__cpluplus -O ' ];

ICR_OBJECTS =	[ 'icr.cpp ', '../lib/libicr.a ', '../tools/obj-1.0/obj/.libs/libobj.so ', ...
                  '../tools/Qhull-2011.2/lib/libqhull6.so ', '../tools/Qhull-2011.2/lib/libqhullcpp.a '];

eval( [ 'mex CXX=g++ -output icr ', CPPFLAGS, ICR_OBJECTS ] );

disp( [ 'icr.', eval('mexext'), ' successfully created!'] );

%mex CXX=g++ -output icr icr.cpp -I../tools/eigen-eigen-3.0.2/ -I../tools/Qhull-2011.2/src/libqhullcpp 