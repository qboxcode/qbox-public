#!/bin/bash
# qbox_reduce.sh: Remove wave functions from a restart file
# The original file is modified and contains the <atomset> element only
# use: qbox_reduce.sh file.xml [file.xml ..]
for f in ${*}
do
  tmpfile=qbox_reduce$$
  nlines=$(grep /atomset -m 1 -n $f | cut -f1 -d: - )
  head -$nlines $f > $tmpfile
  echo "</fpmd:sample>" >> $tmpfile
  mv $tmpfile $f
done
