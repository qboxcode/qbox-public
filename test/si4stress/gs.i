# Si4 ground state
set cell 16 0 0  0 16 0  0 0 16
set ref_cell 18 0 0  0 18 0  0 0 18
species silicon Si_VBC_LDA-1.0.xml
atom Si1 silicon  3.700 -0.100  0.300
atom Si2 silicon -0.100  2.800 -0.200
atom Si3 silicon -3.700  0.100  0.300
atom Si4 silicon -0.100 -2.800 -0.200
strain 0.02 0.04 0.06 0.03 0.05 0.07
params.i
set wf_dyn PSDA
set stress ON
set debug STRESS ON
randomize_wf
set scf_tol 1.e-9
run 0 300
save gs.xml
