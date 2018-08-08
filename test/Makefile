include Makefile.inc
testdirs:= h2ogs h2ocg si2gs cgcell si4gs si4cg si4mlwf si4bomd si4cpmd 
           efield o2gs si2stress si4stress
all:
	for d in $(testdirs); do (cd $$d; $(MAKE)); done
clean:
	for d in $(testdirs); do (cd $$d; $(MAKE) clean); done
