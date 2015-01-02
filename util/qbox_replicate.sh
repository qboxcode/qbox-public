#!/bin/bash
#
# qbox_replicate: replicate the unit cell in the a0,a1,a2 directions.
# use: qbox_replicate cell.sys n0 n1 n2 > newcell.sys
#
if (( $# != 4 ))
  then echo "use: qbox_replicate cell.sys n0 n1 n2 > newcell.sys"
  exit
fi
gawk -v n0=$2 -v n1=$3 -v n2=$4 \
    ' / cell/  {a0x=$3;a0y=$4;a0z=$5; \
               a1x=$6;a1y=$7;a1z=$8; \
               a2x=$9;a2y=$10;a2z=$11; \
      print "set cell ", \
        n0*$3,n0*$4,n0*$5, n1*$6,n1*$7,n1*$8, n2*$9,n2*$10,n2*$11} \
      /ref_cell/  { \
      print "set ref_cell ", \
        n0*$3,n0*$4,n0*$5, n1*$6,n1*$7,n1*$8, n2*$9,n2*$10,n2*$11} \
      /species/  {print} \
      /atom/ {x=$4 - (n0-1)*a0x/2 - (n1-1)*a1x/2 - (n2-1)*a2x/2; \
              y=$5 - (n0-1)*a0y/2 - (n1-1)*a1y/2 - (n2-1)*a2y/2; \
              z=$6 - (n0-1)*a0z/2 - (n1-1)*a1z/2 - (n2-1)*a2z/2; \
              for ( i=0; i<n0; i++ )
              for ( j=0; j<n1; j++ )
              for ( k=0; k<n2; k++ )
                printf("atom %s_%d%d%d %s %12.6f %12.6f %12.6f\n", \
                       $2,i,j,k,$3, \
                  x+i*a0x+j*a1x+k*a2x, \
                  y+i*a0y+j*a1y+k*a2y, \
                  z+i*a0z+j*a1z+k*a2z) \
             }' $1
