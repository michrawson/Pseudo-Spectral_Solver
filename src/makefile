#-std=f95
FLGS=-I. -I/home/mr4209/github/Pseudo-Spectral_Solver/fftw-3.3.4/api/ -fPIC
#FLGS+=-Ofast -march=native
FLGS+=-O3 -march=native
FLGS+=-g -fbounds-check -fbacktrace -Wall -Wextra -Wconversion -pedantic 

SRC=$(wildcard *.f95)
OBJS=$(SRC:.f95=.o)

LIB_SRC=$(wildcard lib/*.f95)
LIB_OBJS=$(LIB_SRC:.f95=.o)

main.tsk : $(LIB_OBJS) $(OBJS)
	( cd ./Shidong_Poisson; $(MAKE) libnufft.a; $(MAKE) poisson2df.o )
	gfortran $(FLGS) \
                Shidong_Poisson/poisson2df.o $(LIB_OBJS) $(OBJS) Shidong_Poisson/*.a -o main.tsk \
        		-L/home/mr4209/github/Pseudo-Spectral_Solver/fftw-3.3.4/.libs/ \
                -lfftw3

fmain : Shidong_Poisson/poisson2df.o
	f2py -c --debug-capi --fcompiler=gfortran \
                -I/home/mr4209/github/Pseudo-Spectral_Solver/fftw-3.3.4/api/ \
				-L../fftw-3.3.4/.libs \
                -lfftw3 \
                -m fmain \
                lib/*.f95 Shidong_Poisson/poisson2df.o Shidong_Poisson/*.a 

%.o : %.f95
	gfortran -c $(FLGS) $< -o $@

%.o : %.f
	gfortran -c $(FLGS) $< -o $@

clean:
	( cd ./Shidong_Poisson; $(MAKE) clean )
	( cd lib; rm -rf *.mod *.o *.so *.tsk *.tsk.dSYM )
	rm -rf *.mod *.o *.so *.tsk *.tsk.dSYM *.so.dSYM

