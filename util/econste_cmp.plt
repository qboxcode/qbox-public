#!/bin/bash
# econste_cmp.plt: compare <econst> and <etotal> in two MD simulations
# use: econste_cmp.plt mdrun1.r mdrun2.r
if [ $1 == "-range" ]
then
  range=$2
  shift 2
fi
gnuplot -persist <<EOF
p $range "<grep -h econst $1" u 2 w l, "<grep -h '<etotal>' $1" u 2 w l, \
  "<grep -h econst $2" u 2 w l, "<grep -h '<etotal>' $2" u 2 w l

EOF
