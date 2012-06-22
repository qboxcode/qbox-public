# Si4 geometry optimization
load ../si4gs/test.xml
set wf_dyn PSDA
set ecutprec 2
set atoms_dyn CG
run 50 5
save test.xml
