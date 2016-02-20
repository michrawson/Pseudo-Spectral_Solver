#-std=f95
FLGS=-g -I. -fbounds-check -fbacktrace \
		 -Wall -Wextra -Wconversion -pedantic

main.tsk :
	gfortran $(FLGS) \
				-L/Users/mich/github/Pseudo-Spectral_Solver/fftpack5.1/lib \
				-L/Users/mich/github/Pseudo-Spectral_Solver/fftw-3.3.4/.libs \
                -lfftpack -lfftw3 \
                lib/*.f95 *.f95  -o main.tsk
	f2py --fcompiler=gfortran --debug-capi \
				-L/Users/mich/github/Pseudo-Spectral_Solver/fftpack5.1/lib \
				-L/Users/mich/github/Pseudo-Spectral_Solver/fftw-3.3.4/.libs \
                -lfftpack -lfftw3 \
                -c lib/*.f95 *.f95 \
                -m fmain

clean:
	rm -rf *.mod *.o *.so *.tsk *.tsk.dSYM
