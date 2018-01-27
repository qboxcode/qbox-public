load gs.xml
params.i
set wf_dyn PSDA
set stress ON
set debug STRESS
run 0
set wf_dyn LOCKED
strain 0 0 0  0.001 0.002 -0.001
run 0
strain -inverse 0 0 0  0.001 0.002 -0.001
strain -inverse 0 0 0  0.001 0.002 -0.001
run 0
