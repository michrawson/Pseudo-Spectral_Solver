
FLGS=-g -std=f2008 -I. -fbounds-check -fbacktrace \
		 -Wall -Wextra -Wconversion -pedantic

main.tsk :
	gfortran $(FLGS) \
				-L/Users/mich/github/Pseudo-Spectral_Solver/fftpack5.1/lib \
                -lfftpack \
                *.f08  -o main.tsk

clean:
	rm -f *.tsk
