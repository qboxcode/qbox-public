set cell        8 0 0  0 8 0  0 0 8
species carbon C_HSCV_PBE-1.0.xml
species hydrogen H_HSCV_PBE-1.0.xml
atom C    carbon       0.00000000   0.00000000   0.00000000
atom H1   hydrogen     1.20000000   1.20000000   1.20000000
atom H2   hydrogen     1.20000000  -1.20000000  -1.20000000
atom H3   hydrogen    -1.20000000   1.20000000  -1.20000000
atom H4   hydrogen    -1.20000000  -1.20000000   1.20000000
set ecut 35
set xc PBE
randomize_wf
set wf_dyn PSDA
set ecutprec 5
run 1 200
