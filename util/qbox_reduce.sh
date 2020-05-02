#!/bin/bash
# reduce qbox restart file to atomset file
# use: qbox_reduce.sh file.xml [file.xml ..]
for f in ${*}
do
  name=${f%.xml}
  atomset_name=${name}_atomset
  echo $name.xml "->" $atomset_name.xml
  get_atomset $f > $atomset_name.xml
  rm $f
done
