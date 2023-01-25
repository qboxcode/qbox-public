#!/bin/bash
if [ $1 == "-range" ]
then
  range=$2
  shift 2
fi
gnuplot -persist <<EOF
set xrange "step"
set ylabel "energy (Ha)"
p $range "<grep -h econst $*" u 2 w l, "<grep -h '<etotal>' $*" u 2 w l
EOF
