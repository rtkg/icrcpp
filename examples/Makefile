#############################################################################
# Makefile for building the example demonstrating libicr.a
# Robert Krug 2011-08-10
#############################################################################

CC = g++ -g -pthread
#CFLAGS = -pipe -Wall -W  -Wextra -Wwrite-strings -Wno-sign-conversion -O3 -g -fPIC -msse2 -std=gnu++0x 
CFLAGS = -O3 -g -fPIC -std=gnu++0x 

LIBS  = -L ../lib -licr                                   #include libicr library
LIBS += -Wl,-rpath,../tools/obj-1.0/build/lib -L ../tools/obj-1.0/build/lib -lobj               #include libobj library
LIBS += -Wl,-rpath,../tools/Qhull-2011.2/lib -L ../tools/Qhull-2011.2/lib -lqhullcpp -lqhull6  #include qhull libs LINK ORDER MATTERS!!!
LIBS += -Wl,-rpath,/home/rkg/Data/Programs/gurobi562/linux64/lib -L/home/rkg/Data/Programs/gurobi562/linux64/lib/ -lgurobi_c++ -lgurobi56 

FN=experiments_ICR_tilt
#FN=experiments_ICR_shift
#FN=example

IDIR += -I ../tools/eigen-eigen-3.0.2            #include eigen headers
IDIR += -I ../tools/Qhull-2011.2/src             #include qhull headers
IDIR += -I ../tools/obj-1.0/build/include        #include libobj headers
IDIR += -I /home/rkg/Data/Programs/gurobi562/linux64/include #Gurobi headers

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
