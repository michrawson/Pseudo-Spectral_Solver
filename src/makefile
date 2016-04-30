#-std=f95
FLGS=-I. -freal-4-real-8 -I../fftw-3.3.4/api/ -fPIC
FLGS+=-Ofast -march=native
#FLGS+=-O3 -march=native
#FLGS+=-g -fbounds-check -fbacktrace -Wall -Wextra  -pedantic  -fcheck=all -Warray-temporaries -fimplicit-none -ffree-line-length-0 -ffpe-trap=zero,overflow -finit-real=nan      -Wconversion     

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


main.tsk : $(aOBJS) $(bOBJS) $(cOBJS)
	( cd ./Shidong_Poisson; $(MAKE) libnufft.a; $(MAKE) poisson2df.o )
	gfortran $(FLGS) \
                Shidong_Poisson/poisson2df.o $(aOBJS) $(bOBJS) $(cOBJS) Shidong_Poisson/*.a -o main.tsk \
        		-L../fftw-3.3.4/.libs/ \
                -lfftw3

fmain : $(aOBJS) $(bOBJS) $(cOBJS) Shidong_Poisson/poisson2df.o
	rm -f lib/variable_coeff_wave_eq-pseudo-RK.o
	f2py --fcompiler=gfortran \
                -I../fftw-3.3.4/api/ \
				-L../fftw-3.3.4/.libs \
                -lfftw3 \
                -c c_fmain.f95 $(aOBJS) $(bOBJS) Shidong_Poisson/poisson2df.o Shidong_Poisson/*.a \
                -m fmain
# --debug-capi

%.o : %.f95
	gfortran -c $(FLGS) $< -o $@

%.o : %.f
	gfortran -c $(FLGS) $< -o $@

clean:
	( cd ./Shidong_Poisson; $(MAKE) clean )
	( cd lib; rm -rf *.mod *.o *.so *.tsk *.tsk.dSYM )
	rm -rf *.mod *.o *.so *.tsk *.tsk.dSYM *.so.dSYM

