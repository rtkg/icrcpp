IFLAGS  = [ '-I. -I../tools/eigen-eigen-3.0.2/ ' ];

CPPFLAGS  = [ IFLAGS, '-D__cpluplus -O ' ];

ICR_OBJECTS =	[ 'icr.cpp '];

eval( [ 'mex CXX=g++ -output icr ', CPPFLAGS, ICR_OBJECTS ] );

disp( [ 'icr.', eval('mexext'), ' successfully created!'] );