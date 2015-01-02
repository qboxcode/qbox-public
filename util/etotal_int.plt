#!/bin/bash
# etotal_int.plt: plot <etotal_int> in one or more simulations
# use: etotal_int.plt run1.r  [run2.r ...]
gnuplot -persist <<EOF
plot "<grep -h etotal $*" u 2 w l
EOF
