#!/bin/bash
# etotal.plt: plot <etotal> in one or more simulations
# use: etotal.plt run1.r [run2.r ...]
gnuplot -persist <<EOF
plot "<grep -h '<etotal>' $*" u 2 w l
EOF
