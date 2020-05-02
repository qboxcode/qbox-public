#!/bin/bash
# qbox_translate: translate all atoms
# use: qbox_replicate cell.sys dx dy dz > newcell.sys
#
if (( $# != 4 ))
  then echo "use: qbox_translate cell.sys dx dy dz > newcell.sys"
  exit
fi
awk -v dx=$2 -v dy=$3 -v dz=$4 \
    ' /^#/  {print} \
      / cell/  {print} \
      /ref_cell/  {print}
      /species/  {print} \
      /atom/ {x=$4 + dx; y=$5 + dy; z=$6 + dz; \
              printf("atom %s %s %12.6f %12.6f %12.6f\n", \
                       $2,$3,x,y,z) \
             }' $1
