include Makefile.inc
testdirs:= h2ogs si2gs cgcell h2ocg
all:
	for d in $(testdirs); do (cd $$d; $(MAKE)); done
clean:
	for d in $(testdirs); do (cd $$d; $(MAKE) clean); done
