#-std=f95
FLGS=-g -I. -fbounds-check -fbacktrace \
		 -Wall -Wextra -Wconversion -pedantic -fdefault-real-8 

main.tsk :
	gfortran $(FLGS) \
				-L../fftpack5.1/lib \
				-L../fftw-3.3.4/.libs \
                -lfftpack -lfftw3 \
                lib/*.f95 *.f95  -o main.tsk
	f2py --fcompiler=gfortran --debug-capi \
				-L../fftpack5.1/lib \
				-L../fftw-3.3.4/.libs \
                -lfftpack -lfftw3 \
                -c lib/*.f95 \
                -m fmain
	# *.f95 \

clean:
	rm -rf *.mod *.o *.so *.tsk *.tsk.dSYM
