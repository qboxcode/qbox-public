#!/bin/bash
if [ $1 == "-range" ]
then
  range=$2
  shift 2
fi
gnuplot -persist <<EOF
set grid
plot $range "<grep -h '<unit_cell_volume>' $*" u 2 w l
EOF
