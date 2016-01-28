FC = gfortran
FLGS = -g -pg -std=f2008 -I. -fbounds-check
#DEPS = camx.prm
#OBJ = helloworld.o

#linker macro
# %.o: %.F08 $(DEPS)

#build targets

finite_differences.tsk:
	        find . -maxdepth 1 -type f -name '*.[fF]??' | xargs $(FC) -o $@ $(FLGS)
			find . -maxdepth 1 -type f -name '*.[fF]??' | xargs -I {} f2py -c {} -m finite_differences 

#	        find . -maxdepth 1 -type f -name '*.[fF]??' | tr '\n' ' ' | $(FC) -o $@ $< $(FLGS)
#	        find . -maxdepth 1 -type f -name '*.[fF]??' | tr '\n' ' ' | $(FC) *.{f,F}?? -o $@ $< $(FLGS)
#	        $(FC)    *.o -o $@ $^ $(FLGS)

clean:
	        rm -fr *.so *.o *.tsk

