#!/bin/bash
# econste.plt: plot <econst> and <etotal> in one or more MD simulations
# use: econste.plt mdrun1.r [mdrun2.r ...]
gnuplot -persist <<EOF
p "<grep -h econst $*" u 2 w l, "<grep -h '<etotal>' $*" u 2 w l
EOF
