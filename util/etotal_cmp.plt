#!/bin/bash
if [ $1 == "-range" ]
then
  range=$2
  shift 2
fi
gnuplot -persist <<EOF
set xlabel "step"
set ylabel "energy (Ha)"
plot $range "<grep '<etotal>' $1" u 2 w l, "<grep '<etotal>' $2" u 2 w l
EOF
