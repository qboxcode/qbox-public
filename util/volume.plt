#!/bin/bash
# volume.plt: plot the unit cell volume during one or more MD simulations
# use: volume.plt mdrun1.r [mdrun2.r ...]
gnuplot -persist <<EOF
set grid
plot "<grep -h '<unit_cell_volume>' $*" u 2 w l
EOF
