set cell 16 0 0  0 16 0  0 0 16
species oxygen O_HSCV_PBE-1.0.xml
species hydrogen H_HSCV_PBE-1.0.xml
atom O oxygen 0 0 0
atom H1 hydrogen  1.128  1.444  0.000
atom H2 hydrogen  1.128 -1.444  0.000

distance O H1
distance O H2
angle H1 O H2

set ecut 70
set wf_dyn PSDA
set ecutprec 5
randomize_wf
run 0 100
save test.xml
