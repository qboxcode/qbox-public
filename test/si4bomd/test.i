# Si4 BO dynamics
load ../si4gs/test.xml
set wf_dyn PSDA
set ecutprec 2
set atoms_dyn MD
set dt 40
run 50 5
