#-std=f95
FFLAGS=-I. -freal-4-real-8 -I../fftw-3.3.4/api/ -fPIC
FFLAGS+=-Ofast -march=native
#FFLAGS+=-O0 -march=native -g -fbounds-check -fbacktrace -Wall -Wextra  -pedantic  -fcheck=all -Warray-temporaries -fimplicit-none -ffree-line-length-0 -ffpe-trap=zero,overflow -finit-real=nan      -Wconversion     

aSRC=$(wildcard a_*.f95)
aOBJS=$(aSRC:.f95=.o)

bSRC=$(wildcard b_*.f95)
bOBJS=$(bSRC:.f95=.o)

#fmain : Shidong_Poisson/poisson2df.o
#	f2py -c --debug-capi --fcompiler=gfortran \
#                -I/home/mr4209/github/Pseudo-Spectral_Solver/fftw-3.3.4/api/ \
#				-L../fftw-3.3.4/.libs \
 #               -lfftw3 \
 #               -m fmain \
 #               lib/*.f95 Shidong_Poisson/poisson2df.o Shidong_Poisson/*.a 

cSRC=$(wildcard c_*.f95)
cOBJS=$(cSRC:.f95=.o)

.PHONY: main.tsk fmain.so

main.tsk : $(aOBJS) $(bOBJS) $(cOBJS)
	( cd ./Shidong_Poisson; $(MAKE) libnufft.a; $(MAKE) poisson2df.o )
	gfortran $(FFLAGS) \
                Shidong_Poisson/poisson2df.o $(aOBJS) $(bOBJS) $(cOBJS) Shidong_Poisson/*.a -o main.tsk \
        		-L../fftw-3.3.4/.libs/ \
                -lfftw3

fmain.so : $(aOBJS) $(bOBJS)
	( cd ./Shidong_Poisson; $(MAKE) libnufft.a; $(MAKE) poisson2df.o )
	f2py --fcompiler=gfortran \
                -I../fftw-3.3.4/api/ \
				-L../fftw-3.3.4/.libs \
                -lfftw3 \
                -c c_fmain.f95 $(aOBJS) $(bOBJS) Shidong_Poisson/poisson2df.o \
                Shidong_Poisson/libnufft.a \
                -m fmain
# --debug-capi

%.o : %.f95
	gfortran -c $(FFLAGS) $< -o $@

%.o : %.f
	gfortran -c $(FFLAGS) $< -o $@

clean:
	( cd ./Shidong_Poisson; $(MAKE) clean )
	( cd lib; rm -rf *.mod *.o *.so *.tsk *.tsk.dSYM )
	rm -rf *.mod *.o *.so *.tsk *.tsk.dSYM *.so.dSYM

