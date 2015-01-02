#!/bin/bash
# econste_cmp.plt: compare <econst> and <etotal> in two MD simulations
# use: econste_cmp.plt mdrun1.r mdrun2.r
gnuplot -persist <<EOF
p "<grep -h econst $1" u 2 w l, "<grep -h '<etotal>' $1" u 2 w l, \
  "<grep -h econst $2" u 2 w l, "<grep -h '<etotal>' $2" u 2 w l

EOF
