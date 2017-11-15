#!/bin/bash
# cell.plt: plot cell parameters during an MD simulation
if [ $1 == "-range" ]
then
  range=$2
  shift 2
fi

gnuplot -persist <<EOF
plot $range "<grep -h -A3 '<unit_cell' $* | grep a=" u 2 w l, \
     "<grep -h -A3 '<unit_cell' $* | grep a=" u 3 w l, \
     "<grep -h -A3 '<unit_cell' $* | grep a=" u 4 w l, \
     "<grep -h -A3 '<unit_cell' $* | grep b=" u 2 w l, \
     "<grep -h -A3 '<unit_cell' $* | grep b=" u 3 w l, \
     "<grep -h -A3 '<unit_cell' $* | grep b=" u 4 w l, \
     "<grep -h -A3 '<unit_cell' $* | grep c=" u 2 w l, \
     "<grep -h -A3 '<unit_cell' $* | grep c=" u 3 w l, \
     "<grep -h -A3 '<unit_cell' $* | grep c=" u 4 w l, 0
fit ax "<grep -h -A3 '<unit_cell' $* | grep a=" u 0:2 via ax
fit ay "<grep -h -A3 '<unit_cell' $* | grep a=" u 0:3 via ay
fit az "<grep -h -A3 '<unit_cell' $* | grep a=" u 0:4 via az
fit bx "<grep -h -A3 '<unit_cell' $* | grep b=" u 0:2 via bx
fit by "<grep -h -A3 '<unit_cell' $* | grep b=" u 0:3 via by
fit bz "<grep -h -A3 '<unit_cell' $* | grep b=" u 0:4 via bz
fit cx "<grep -h -A3 '<unit_cell' $* | grep c=" u 0:2 via cx
fit cy "<grep -h -A3 '<unit_cell' $* | grep c=" u 0:3 via cy
fit cz "<grep -h -A3 '<unit_cell' $* | grep c=" u 0:4 via cz
print "avg a: ", ax, ay, az
print "avg b: ", bx, by, bz
print "avg c: ", cx, cy, cz
EOF
