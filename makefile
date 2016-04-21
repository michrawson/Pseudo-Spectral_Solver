
all: libfftw libfftpack pss

libfftw :
	( cd ./fftw-3.3.4; $(MAKE) )

libfftpack :
	( cd ./fftpack5.1; $(MAKE) )

pss :
	( cd ./src; $(MAKE) )

clean:
	( cd ./fftpack5.1; $(MAKE) clean; cd ../fftw-3.3.4; $(MAKE) clean; cd ../src; $(MAKE) clean )
