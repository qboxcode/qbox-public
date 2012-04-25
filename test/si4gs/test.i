# Si4 ground state
set cell 16 0 0  0 16 0  0 0 16
species silicon Si_VBC_LDA-1.0.xml
atom Si1 silicon  3.700  0.000  0.000
atom Si2 silicon  0.000  2.200  0.000
atom Si3 silicon -3.700  0.000  0.000
atom Si4 silicon  0.000 -2.200  0.000

set ecut 6
set wf_dyn PSDA
set ecutprec 2

randomize_wf
run 0 200
save test.xml
