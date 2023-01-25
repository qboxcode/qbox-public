#!/bin/bash
if [ $1 == "-range" ]
then
  range=$2
  shift 2
fi
gnuplot -persist <<EOF
set xrange "step"
set ylabel "energy (Ha)"
plot $range "<grep -h '<etotal>' $*" u 2 w l
EOF
