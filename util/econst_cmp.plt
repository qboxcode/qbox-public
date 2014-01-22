#!/bin/bash
# econst_cmp.plt: compare <econst> in two MD simulations
# use: econst_cmp.plt mdrun1.r mdrun2.r
gnuplot -persist <<EOF
plot "<grep '<econst>' $1" u 2 w l, "<grep '<econst>' $2" u 2 w l
EOF
