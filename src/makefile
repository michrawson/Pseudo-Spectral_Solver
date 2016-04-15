#-std=f95
FLGS=-O2 -I. -freal-4-real-8
FLGS+=-g -fbounds-check -fbacktrace -Wall -Wextra -Wconversion -pedantic -fcheck=all

SRC=$(wildcard *.f95)
OBJS=$(SRC:.f95=.o)

MOD_SRC=$(wildcard lib/*.f95)
MODS=$(MOD_SRC:.f95=.mod)

main.tsk : $(MODS) $(OBJS) 
	( cd ./Shidong_Poisson; $(MAKE) nufftlib; $(MAKE) clean; $(MAKE) poisson2df )
	gfortran $(FLGS) \
				-L../fftpack5.1/lib \
				-L../fftw-3.3.4/.libs \
                -lfftpack -lfftw3 \
                Shidong_Poisson/*.o $(MODS) $(OBJS) Shidong_Poisson/*.a -o main.tsk

fmain : $(SRC)
	f2py --fcompiler=gfortran --debug-capi \
				-L../fftpack5.1/lib \
				-L../fftw-3.3.4/.libs \
                -lfftpack -lfftw3 \
                -c lib/*.f95 Shidong_Poisson/*.o Shidong_Poisson/*.a \
                -m fmain

%.mod : %.f95
	gfortran -c $(FLGS) $< -o $@

%.o : %.f95
	gfortran -c $(FLGS) $< -o $@

clean:
	( cd ./Shidong_Poisson; $(MAKE) clean )
	rm -rf lib/*.mod lib/*.o *.mod *.o *.so *.tsk *.tsk.dSYM
