#!/bin/bash
# qbox_rotate.sh rotate sample counterclockwise by $1 degrees around
# x, y, or z axis
# use: qbox_rotate.sh {x|y|z} angle file.i
if [ $# != 3 ]
then
  echo 'use: qbox_rotate.sh {x|y|z} angle file.i'
  exit
fi

axis=$1
theta=$2
c=$(bc -l <<< "c($theta * 4 * a(1.0) / 180)")
s=$(bc -l <<< "s($theta * 4 * a(1.0) / 180)")

if [ $axis == 'x' ]
then
  awk -v c=$c -v s=$s \
  '/^#/ {print $0} \
  /species/ {print $0} \
  /set cell/ {ax=$3; ay=$4; az=$5; \
              bx=$6; by=$7; bz=$8; \
              cx=$9; cy=$10; cz=$11; \
    printf("%s %12.6f %12.6f %12.6f  \
                  %12.6f %12.6f %12.6f  \
                  %12.6f %12.6f %12.6f\n", \
            "set cell ", \
            ax,c*ay-s*az,s*ay+c*az, \
            bx,c*by-s*bz,s*by+c*bz, \
            cx,c*cy-s*cz,s*cy+c*cz)} \
  /atom/ {x=$4; y=$5; z=$6; printf("%s %s %s %12.6f %12.6f %12.6f\n", \
  $1,$2,$3,x,c*y-s*z,s*y+c*z)}' $3
elif [ $axis == 'y' ]
then
  awk -v c=$c -v s=$s \
  '/^#/ {print $0} \
  /species/ {print $0} \
  /set cell/ {ax=$3; ay=$4; az=$5; \
              bx=$6; by=$7; bz=$8; \
              cx=$9; cy=$10; cz=$11; \
    printf("%s %12.6f %12.6f %12.6f  \
                  %12.6f %12.6f %12.6f  \
                  %12.6f %12.6f %12.6f\n", \
            "set cell ", \
            c*ax+s*az,ay,-s*ax+c*az, \
            c*bx+s*bz,by,-s*bx+c*bz, \
            c*cx+s*cz,cy,-s*cx+c*cz)} \
  /atom/ {x=$4; y=$5; z=$6; printf("%s %s %s %12.6f %12.6f %12.6f\n", \
  $1,$2,$3,c*x+s*z,y,-s*x+c*z)}' $3
elif [ $axis == 'z' ]
then
  awk -v c=$c -v s=$s \
  '/^#/ {print $0} \
  /species/ {print $0} \
  /set cell/ {ax=$3; ay=$4; az=$5; \
              bx=$6; by=$7; bz=$8; \
              cx=$9; cy=$10; cz=$11; \
    printf("%s %12.6f %12.6f %12.6f  \
                  %12.6f %12.6f %12.6f  \
                  %12.6f %12.6f %12.6f\n", \
            "set cell ", \
            c*ax-s*ay,s*ax+c*ay,az, \
            c*bx-s*by,s*bx+c*by,bz, \
            c*cx-s*cy,s*cx+c*cy,cz)} \
  /atom/ {x=$4; y=$5; z=$6; printf("%s %s %s %12.6f %12.6f %12.6f\n", \
  $1,$2,$3,c*x-s*y,s*x+c*y,z)}' $3
else
  echo "invalid axis name: " $axis
  echo 'use: qbox_rotate.sh {x|y|z} angle file.i'
fi
