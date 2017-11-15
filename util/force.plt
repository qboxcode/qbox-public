#!/bin/bash
if [ $1 == "-range" ]
then
  range=$2
  shift 2
fi
gnuplot -persist <<EOF
p $range "<grep -h force $*" u 2, "" u 3, "" u 4, 0
EOF
