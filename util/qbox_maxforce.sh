#!/bin/bash
# get the largest force component in a file of "<force> fx fy fz </force>"
# use: qbox_maxforce nat file.r
grep '<force>' $2 | tail -$1 | \
awk '{
  if ( mx*mx < $2*$2 ) mx = $2;
  if ( my*my < $3*$3 ) my = $3;
  if ( mz*mz < $4*$4 ) mz = $4;
} END {printf("%9.2e %9.2e %9.2e\n", mx, my, mz)}' - 
