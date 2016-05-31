#!/bin/bash
if [ $1 == "-range" ]
then
  range=$2
  shift 2
fi
gnuplot -persist <<EOF
plot $range "<grep -h eigenvalue_sum $*" u 2 w l
EOF
