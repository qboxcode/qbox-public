include ../Makefile.inc
all: test.r
clean:
	rm -f gs.r gs.xml gs_diag.r gs_offdiag.r test.r

gs.r: params.i
gs_diag.r: params.i
gs_offdiag.r: params.i

gs.r: gs.i
	$(QBOX_CMD) gs.i > gs.r

gs_diag.r: gs_diag.i gs.r
	$(QBOX_CMD) gs_diag.i > gs_diag.r

gs_offdiag.r: gs_offdiag.i gs.r
	$(QBOX_CMD) gs_offdiag.i > gs_offdiag.r

test.r: gs_diag.r gs_offdiag.r
	./check_stress.sh gs_diag.r > test.r
	./check_stress.sh gs_offdiag.r >> test.r
