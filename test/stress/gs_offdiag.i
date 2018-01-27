load gs.xml
params.i
set wf_dyn PSDA
set stress ON
set debug STRESS
run 0
set wf_dyn LOCKED
strain 0 0 0  0.0005 0.000  0.000
run 0
strain -inverse 0 0 0  0.0005 0.000  0.000
strain -inverse 0 0 0  0.0005 0.000  0.000
run 0
