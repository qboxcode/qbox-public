#!/bin/bash
# etotal_cmp.plt: compare <etotal> in two MD simulations
# use: etotal_cmp.plt mdrun1.r mdrun2.r
gnuplot -persist <<EOF
set grid
plot "<grep '<etotal>' $1" u 2 w l, "<grep '<etotal>' $2" u 2 w l
EOF
