include Makefile.inc
testdirs:= h2ogs h2ocg h2oresp \
           si2gs si2stress cgcell \
           si4gs si4cg si4mlwf si4bomd si4cpmd  \
           efield o2gs si4stress siatom
all:
	for d in $(testdirs); do $(MAKE) -C $$d -w; done
clean:
	for d in $(testdirs); do $(MAKE) -C $$d -w clean; done
