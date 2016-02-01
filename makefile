
include make.inc

all: libfftpack testfftpack

libfftpack:
	mkdir -p ./objs
	( cd ./src; $(MAKE) clean; $(MAKE) )

clean:
	( cd ./src; $(MAKE) clean )
