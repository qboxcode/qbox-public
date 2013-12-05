#!/bin/bash
# force.plt: plot all components of ionic forces in one or more simulations
# use: force.plt mdrun1.r [mdrun2.r ...]
gnuplot -persist <<EOF
p "<grep -h force $*" u 2, "" u 3, "" u 4, 0
EOF
