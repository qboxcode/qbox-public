#!/bin/bash
# compute Egap = E(n+1) - E(n)
# use: egap.sh n file.r
declare -i nocc=$1
shift
declare -i nl=nocc/5+1
#echo "nl=" $nl
declare -i nfrac=nocc-5*\(nocc/5\)
#echo "nfrac=" $nfrac

grep -h -A$nl '<eigenvalues ' ${*}| \
 awk -v nl=$nl -v nfrac=$nfrac \
  ' NR%(nl+2)==nl {e[0] = $5;} \
    NR%(nl+2)==(nl+1) \
    { e[1]=$1; e[2]=$2; e[3]=$3; e[4]=$4; e[5]=$5; \
      e_homo = e[nfrac]; e_lumo = e[nfrac+1]; \
      print "E_HOMO=",e_homo, "E_LUMO=",e_lumo, "Eg=",e_lumo-e_homo}' -
