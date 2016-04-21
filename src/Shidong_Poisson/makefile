all: nufftlib int2

PROJECT=int2         # 2d poisson

RM      = /bin/rm

CC	= gcc

#stheno
#FC      = ifort

#DEBUG	= -quad

#FFLAGS  =  -O3 -c

#LIBDIR	= /opt/fftw3/3.1.2/intel/lib

#LIBS 	= -L$(LIBDIR) -lfftw3 -lm

#LINKER 	= ifort

#access
FC      = gfortran

#DEBUG	= -quad

FFLAGS  =  -O3 -c 

# -g -fbacktrace -fbounds-check -fcheck=all

LIBDIR=../../fftw-3.3.4/.libs

LIBS 	= -L$(LIBDIR) -lfftw3 -lm
# -lfftw3f

# laptop
#LIBDIR=

#LIBS 	= -lfftw3f -lfftw3 

LINKER	= gfortran

AR	= ar

ARFLAGS	= -r

RANLIB	= ranlib

NUFFTLIB = libnufft.a

TARGET=int2

vpath %.f .:./common
vpath %.f90 .:./common

AROBJS 	=dfftpack.o next235.o nufft2df90.o\
	second-r8.o legeexps.o prini.o hkrand.o dlaran.o\
	eone.o erfc.o besei0.o besei1.o calcei.o calci0.o calci1.o

libnufft.a :  $(AROBJS)
	rm -f $(NUFFTLIB)
	$(AR) $(ARFLAGS) $(NUFFTLIB) $(AROBJS) 
	$(RANLIB) $(NUFFTLIB)

$(PROJECT): testpoisson2df.o poisson2df.o
	rm -f $(PROJECT)
	$(LINKER) $(DEBUG) -o $(PROJECT) testpoisson2df.o poisson2df.o libnufft.a $(LIBS)
#	./$(PROJECT)

poisson2df: poisson2df.o
	:

#$(PROJECT): testpoisson2d.o poisson2d.o
#	rm -f $(PROJECT)
#	$(LINKER) $(DEBUG) -o $(PROJECT) testpoisson2d.o poisson2d.o libnufft.a $(LIBS)
#	./$(PROJECT)

# SUFFIX RULE
.f.o:
	$(FC) $(DEBUG) $(FFLAGS) $<
.f90.o:
	$(FC) $(DEBUG) $(FFLAGS) $<

.SUFFIXES: .o .f .f90

clean:
	$(RM) -f *.o *~ \#* fort.* *.exe $(PROJECT) 

