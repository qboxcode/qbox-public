load ../h2ogs/test.xml
set wf_dyn PSDA
set ecutprec 5

set atoms_dyn CG
run 20 10
distance O H1
distance O H2
angle H1 O H2
save test.xml
