FC = gfortran
FLGS = -g -pg -std=f2008 -I. -fbounds-check
#DEPS = camx.prm
#OBJ = helloworld.o

#linker macro
# %.o: %.F08 $(DEPS)

#build targets

finite_differences.tsk:
	        $(FC) *.F08 -o $@ $< $(FLGS)
#	        $(FC)    *.o -o $@ $^ $(FLGS)

clean:
	        rm -f *.o *.tsk

