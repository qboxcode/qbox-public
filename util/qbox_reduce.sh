#!/bin/bash
# Reduce a qbox restart file to an atomset file
# The original file is removed
# use: qbox_reduce.sh file.xml [file.xml ..]
for f in ${*}
do
  name=${f%.xml}
  atomset_name=${name}_atomset
  echo $name.xml "->" $atomset_name.xml
  nlines=$(grep /atomset -m 1 -n $f | cut -f1 -d: - )
  head -$nlines $f > $atomset_name.xml
  echo "</fpmd:sample>" >> $atomset_name.xml
  rm $f
done
