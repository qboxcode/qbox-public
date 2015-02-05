set nrowmax 3 
# Si4 ground state
set cell 16 0 0  0 16 0  0 0 16
species silicon Si_VBC_LDA-1.0.xml
atom Si1 silicon  3.700  0.000  1.000
atom Si2 silicon  0.000  2.200  1.000
atom Si3 silicon -3.700  0.000  1.000
atom Si4 silicon  0.000 -2.200  1.000

set ecut 6
set wf_dyn PSDA
set ecutprec 2

randomize_wf
run 0 200

set e_field 0 0  0.001
set polarization BERRY
run 0 100
set polarization MLWF_REF
run 0 100
set polarization MLWF_REF_Q
run 0 100
