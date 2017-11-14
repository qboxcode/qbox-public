#!/bin/bash
# generate an xyz file from a k-point input file
awk '/k-points/ { printf("%d\n\n", $2) } \
     /b0/ { b0x=$3; b0y=$4; b0z=$5 } \
     /b1/ { b1x=$3; b1y=$4; b1z=$5 } \
     /b2/ { b2x=$3; b2y=$4; b2z=$5 } \
     /kpoint add/ { \
        kx=$3*b0x+$4*b1x+$5*b2x; \
        ky=$3*b0y+$4*b1y+$5*b2y; \
        kz=$3*b0z+$4*b1z+$5*b2z; \
        printf("H %10.4f %10.4f %10.4f\n",kx,ky,kz)}' $1
