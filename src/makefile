#-std=f95
FLGS=-g -I. -fbounds-check -fbacktrace \
		 -Wall -Wextra -Wconversion -pedantic

main.tsk :
	gfortran $(FLGS) \
				-L/Users/mich/github/Pseudo-Spectral_Solver/fftpack5.1/lib \
                -lfftpack \
                lib/*.f95 *.f95  -o main.tsk

clean:
	rm -rf *.mod *.o *.tsk *.tsk.dSYM
