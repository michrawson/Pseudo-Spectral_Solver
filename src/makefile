#-std=f95
FLGS=-Ofast -march=native -I. -freal-4-real-8
# FLGS+=-g -fbounds-check -fbacktrace -Wall -Wextra -Wconversion -pedantic -fcheck=all

SRC=$(wildcard *.f95)
OBJS=$(SRC:.f95=.o)

LIB_SRC=$(wildcard lib/*.f95)
LIB_OBJS=$(LIB_SRC:.f95=.o)

main.tsk : $(LIB_OBJS) $(OBJS)
	( cd ./Shidong_Poisson; $(MAKE) libnufft.a; $(MAKE) clean; $(MAKE) poisson2df )
	gfortran $(FLGS) \
				-L../fftpack5.1/lib \
				-L../fftw-3.3.4/.libs \
                -lfftpack -lfftw3 \
                Shidong_Poisson/*.o $(LIB_OBJS) $(OBJS) Shidong_Poisson/*.a -o main.tsk

fmain : Shidong_Poisson/poisson2df.o
	rm -f lib/variable_coeff_wave_eq-pseudo-RK.o
	f2py --fcompiler=gfortran --debug-capi \
				-L../fftpack5.1/lib \
				-L../fftw-3.3.4/.libs \
                -lfftpack -lfftw3 \
                -c lib/*.f95 Shidong_Poisson/poisson2df.o Shidong_Poisson/*.a \
                -m fmain

%.o : %.f95
	gfortran -c $(FLGS) $< -o $@

%.o : %.f
	gfortran -c $(FLGS) $< -o $@

clean:
	( cd ./Shidong_Poisson; $(MAKE) clean )
	( cd lib; rm -rf *.mod *.o *.so *.tsk *.tsk.dSYM )
	rm -rf *.mod *.o *.so *.tsk *.tsk.dSYM *.so.dSYM
