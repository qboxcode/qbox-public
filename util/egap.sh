#!/bin/bash
# use: egap.sh n file.r
#
# compute Egap = E(n+1) - E(n)

declare -i nocc=$1
declare -i nl=nocc/5+1
#echo "nl=" $nl
declare -i nfrac=nocc-5*\(nocc/5\)
#echo "nfrac=" $nfrac

grep -A$nl '<eigenvalues ' $2| \
 awk -v nl=$nl -v nfrac=$nfrac -v nocc=$nocc\
  ' /<eigenvalues / {kpx=$3;kpy=$4;kpz=$5} \
    NR%(nl+2)==nl {e[0] = $5;} \
    NR%(nl+2)==(nl+1) \
    { e[1]=$1; e[2]=$2; e[3]=$3; e[4]=$4; e[5]=$5; \
      e_n = e[nfrac]; e_np1 = e[nfrac+1]; \
      printf("E_%d= %-8.3f E_%d= %-8.3f Eg= %-8.3f  %6.3f %6.3f %6.3f\n",nocc,e_n,nocc+1,e_np1,e_np1-e_n,kpx,kpy,kpz);  }' -
