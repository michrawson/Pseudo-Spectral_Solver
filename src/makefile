#-std=f95
FLGS=-g -I. -fbounds-check -fbacktrace \
        -Wall -Wextra -Wconversion -pedantic  -freal-4-real-8 -fcheck=all

main.tsk :
	gfortran -c $(FLGS) `find . -name '*.f95' | sort ` || \
	gfortran -c $(FLGS) `find . -name '*.f95' | sort -r` || \
	gfortran -c $(FLGS) `find . -name '*.f95' | sort ` || \
	gfortran -c $(FLGS) `find . -name '*.f95' | sort -r` || \
	gfortran -c $(FLGS) lib/*.f95 *.f95
	gfortran $(FLGS) \
				-L../fftpack5.1/lib \
				-L../fftw-3.3.4/.libs \
                -lfftpack -lfftw3 \
                *.o -o main.tsk

fmain :
	f2py --fcompiler=gfortran --debug-capi \
				-L../fftpack5.1/lib \
				-L../fftw-3.3.4/.libs \
                -lfftpack -lfftw3 \
                -c lib/*.f95 \
                -m fmain

clean:
	rm -rf *.mod *.o *.so *.tsk *.tsk.dSYM
