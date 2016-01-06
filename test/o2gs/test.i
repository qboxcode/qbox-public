# Oxygen dimer
# triplet state
set cell 11 0 0  0 11 0  0 0 11
species oxygen O_HSCV_PBE-1.0.xml
atom O1 oxygen  1.16  0  0
atom O2 oxygen -1.16  0  0
set nspin 2
set delta_spin 1
set nempty 4
set ecut 70
set wf_dyn JD
randomize_wf
run 0 40 5 
partial_charge O1 1.9
partial_charge -spin 1 O1 1.9
partial_charge -spin 2 O1 1.9
partial_charge O2 1.9
partial_charge -spin 1 O2 1.9
partial_charge -spin 2 O2 1.9
save test.xml
