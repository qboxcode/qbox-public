load gs.xml
params.i
set wf_dyn PSDA
set stress ON
set debug STRESS ON
run 0
set wf_dyn LOCKED
strain 0.0001 0.0002 0.0003  0 0 0
run 0
strain -inverse 0.0001 0.0002 0.0003  0 0 0
strain -inverse 0.0001 0.0002 0.0003  0 0 0
run 0
