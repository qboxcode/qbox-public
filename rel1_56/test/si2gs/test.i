# Si ground state
set cell 5.13 5.13 0 0 5.13 5.13 5.13 0 5.13
set ref_cell 5.5 5.5 0  0 5.5 5.5  5.5 0 5.5
species silicon Si_VBC_LDA-1.0.xml
 atom Si1 silicon 1.2825 1.2825 1.2825
 atom Si2 silicon -1.2825 -1.2825 -1.2825
# Special point in the FCC Brillouin Zone
# A. Baldereschi, Phys. Rev. B7, 5212 (1973)
# 12-point set
kpoint delete 0 0 0
kpoint add  0.45880  0.14765  0.31115   0.0833333333333333
kpoint add -0.16350  0.14765 -0.31115   0.0833333333333333 
kpoint add  0.45880  0.31115  0.14765   0.0833333333333333
kpoint add  0.16350  0.31115 -0.14765   0.0833333333333333
kpoint add  0.31115  0.45880  0.14765   0.0833333333333333
kpoint add -0.31115 -0.16350  0.14765   0.0833333333333333
kpoint add  0.14765  0.45880  0.31115   0.0833333333333333
kpoint add -0.14765  0.16350  0.31115   0.0833333333333333
kpoint add  0.14765  0.31115  0.45880   0.0833333333333333
kpoint add  0.14765 -0.31115 -0.16350   0.0833333333333333
kpoint add  0.31115  0.14765  0.45880   0.0833333333333333
kpoint add  0.31115 -0.14765  0.16350   0.0833333333333333
set stress ON
set debug STRESS
set ecut 15 
set ecuts 10
set wf_dyn PSD
set ecutprec 4
randomize_wf
run 30 10
set wf_diag T
run 0
save test.xml
