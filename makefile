
all: libfftpack pss

libfftpack :
	( cd ./fftpack5.1; $(MAKE) clean; $(MAKE) )

pss :
	( cd ./src; $(MAKE) clean; $(MAKE) )

clean:
	( cd ./fftpack5.1; $(MAKE) clean; cd ../src; $(MAKE) clean )
