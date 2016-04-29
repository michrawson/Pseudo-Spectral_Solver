
all: libfftw pss

libfftw :
	( cd ./fftw-3.3.4; $(MAKE) )

pss :
	( cd ./src; $(MAKE) )

clean:
	( cd fftw-3.3.4; $(MAKE) clean )
	( cd src; $(MAKE) clean )

