#############################################################################
# Makefile for building the example demonstrating libicr.a
# Robert Krug 2011-08-10
#############################################################################

CC = g++ -g -pthread
CFLAGS = -pipe -Wall -W  -Wextra -Wwrite-strings -Wno-sign-conversion -O3 -g -fPIC -msse2 -std=gnu++0x 

LIBS  = -L ../lib -licr                               #include libicr library
LIBS += -L ../tools/obj-1.0/build/lib -lobj           #include libobj library
LIBS += -L ../tools/Qhull-2011.2/lib -lqhull6         #include qhull C library
LIBS += -L ../tools/Qhull-2011.2/lib -lqhullcpp       #include qhull C++ library

FN=example

IDIR += -I ../tools/eigen-eigen-3.0.2            #include eigen headers
IDIR += -I ../tools/Qhull-2011.2/src             #include qhull headers
IDIR += -I ../tools/obj-1.0/build/include        #include libobj headers

E_SDIR = src
E_ODIR = obj

OUT_ODIR = Output

all : ex 
	-@echo SYSTYPE  = $(SYSTYPE)
	
ex : $(E_ODIR)/$(FN).o  
	$(CC)  $< $(IDIR) $(LIBS) -o ex 

$(E_ODIR)/$(FN).o : $(E_SDIR)/$(FN).cpp
	$(CC) -c $< $(CFLAGS) $(IDIR) -o $@

clean :
	-@rm -f $(E_ODIR)/*.o ex*
	-@clear

	
out : 
	-@rm $(OUT_ODIR)/*.dat

###EOF
