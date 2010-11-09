set cell        8 0 0  0 8 0  0 0 8
species carbon http://fpmd.ucdavis.edu/potentials/C/C_HSCV_LDA-1.0.xml
species hydrogen http://fpmd.ucdavis.edu/potentials/H/H_HSCV_LDA-1.0.xml
atom C    carbon       0.00000000   0.00000000   0.00000000
atom H1   hydrogen     1.20000000   1.20000000   1.20000000
atom H2   hydrogen     1.20000000  -1.20000000  -1.20000000
atom H3   hydrogen    -1.20000000   1.20000000  -1.20000000
atom H4   hydrogen    -1.20000000  -1.20000000   1.20000000
set ecut 35
randomize_wf
set wf_dyn PSDA
set ecutprec 5
run 1 200
