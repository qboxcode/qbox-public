# Si atom with custom occupation numbers
set cell 16 0 0  0 16 0  0 0 16
species silicon Si_VBC_LDA-1.0.xml
atom Si silicon  0 0 0
set ecut 6
set nempty 2
set occ 2 0.666666666
set occ 3 0.666666666
set occ 4 0.666666666
randomize_wf
set wf_dyn JD
run 0 30 10
save test.xml
