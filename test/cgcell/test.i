# cgcell.i
# test cell_dyn=CG option

load ../si2gs/test.xml
set stress ON
set debug STRESS
set ecut 15 
set ecuts 10
set wf_dyn PSD
set ecutprec 4

set cell_dyn CG
move Si1 by  0.01 -0.02 -0.005
move Si2 by -0.01  0.02  0.005
strain 0 0 0 0.01 0 0
run 0 30
set cell_dyn CG
run 50 10 
