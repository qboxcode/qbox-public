#!/bin/bash
# etotal_int_cmp.plt: compare <etotal_int> in two simulations
# use: etotal_int_cmp.plt run1.r run2.r
gnuplot -persist <<EOF
plot "<grep '<etotal_int>' $1" u 2 w l, "<grep '<etotal_int>' $2" u 2 w l
EOF
