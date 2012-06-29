include ../Makefile.inc
all: test.r
clean:
	rm -f test.r test.xml

test.r: test.i
	$(QBOX_CMD) $^ > $@
