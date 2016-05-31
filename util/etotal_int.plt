#!/bin/bash
if [ $1 == "-range" ]
then
  range=$2
  shift 2
fi
gnuplot -persist <<EOF
set grid
plot $range "<grep -h etotal_int $*" u 2 w l
EOF
