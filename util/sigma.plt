#!/bin/bash
if [ $1 == "-range" ]
then
  range=$2
  shift 2
fi
gnuplot -persist <<EOF
plot $range "<grep -h sigma_xx $*" u 2 w l, "<grep -h sigma_yy $*" u 2 w l, \
"<grep -h sigma_zz $*" u 2 w l, "<grep -h sigma_xy $*" u 2 w l, \
"<grep -h sigma_yz $*" u 2 w l, "<grep -h sigma_xz $*" u 2 w l, 0
EOF
