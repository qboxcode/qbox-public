#!/bin/bash
if [ $1 == "-range" ]
then
  range=$2
  shift 2
fi
gnuplot -persist <<EOF
fit $range a*x+b "<grep -h temp_ion $*" u 0:2 via a,b
fit $range c "<grep -h temp_ion $*" u 0:2 via c
p $range "<grep -h temp_ion $*" u 2 w l, a*x+b, c

print "Tavg=",c,"   dT/dt=",a
EOF
