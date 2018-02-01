#!/bin/bash
#
# Test the accuracy of the computed stress tensor by comparing it with a
# finite difference calculation of dE=tr(u*sigma)
#
# use: ./check_stress.sh [energy_term sigma_name] file.r
#      If energy_term and sigma_name are not specified,
#      the terms ekin, eps,
# example: ./check_stress.sh ekin sigma_ekin gs.r
#
# The output file file.r must contain the following sequence
# 1) load a sample from a ground state calculation performed
#    using stress=ON, ref_cell, ecuts
# 2) the following commands (with arbitrary strain parameters)
#    set wf_dyn LOCKED
#    strain 0.001 0.001 0.001 0 0 0
#    run 0
#    strain -inverse 0.001 0.001 0.001 0 0 0
#    strain -inverse 0.001 0.001 0.001 0 0 0
#    run 0


function check_stress_term
{
  # check the accuracy of a single term of the energy and stress tensor
  eterm=$1
  etag='<'$eterm'>'
  sigma_name=$2
  file=$3

  echo $eterm $sigma_name

  # read volume
  volume=$(grep '<unit_cell_volume>' $file | awk '{print $2}')
  echo volume: $volume

  # compute energy difference E(+u)-E(-u)
  ep=$(grep '<'$eterm'>' $file | awk 'NR==2 {print $2}')
  em=$(grep '<'$eterm'>' $file | awk 'NR==3 {print $2}')
  echo ep: $ep
  echo em: $em
  de=$(echo "$ep - $em" | bc -l)
  echo de: $de

  # read strain
  u=$(grep '<cmd>' $file | grep strain | head -1 | \
    sed "s/<cmd>//" | sed "s/<\/cmd>//" | sed "s/\[qbox\]//" | sed "s/strain//")
  # echo $u
  strain=($u)
  uxx=$(echo ${strain[0]})
  uyy=$(echo ${strain[1]})
  uzz=$(echo ${strain[2]})
  uxy=$(echo ${strain[3]})
  uyz=$(echo ${strain[4]})
  uxz=$(echo ${strain[5]})
  echo strain: $uxx $uyy $uzz $uxy $uyz $uxz

  # read computed stress tensor
  sxx=$(grep '<'${sigma_name}_xx'>' $file |head -1| awk '{print $2}')
  syy=$(grep '<'${sigma_name}_yy'>' $file |head -1| awk '{print $2}')
  szz=$(grep '<'${sigma_name}_zz'>' $file |head -1| awk '{print $2}')
  sxy=$(grep '<'${sigma_name}_xy'>' $file |head -1| awk '{print $2}')
  syz=$(grep '<'${sigma_name}_yz'>' $file |head -1| awk '{print $2}')
  sxz=$(grep '<'${sigma_name}_xz'>' $file |head -1| awk '{print $2}')
  echo stress: $sxx $syy $szz $sxy $syz $sxz

  # compute energy change from 2 * volume * tr(sigma*u)
  de_comp=$(echo "-2.0 * $volume * ($uxx * $sxx + $uyy * $syy + $uzz * $szz + \
   2.0 * ( $uxy * $sxy + $uyz * $syz + $uxz * $sxz) )" | bc -l)
  printf "E(+u)-E(-u):     %15.6e\n" $de
  printf "2*v*tr(sigma*u): %15.6e\n" $de_comp

  # compute absolute error
  abs_err=$(echo "$de - $de_comp" | bc -l)
  printf "abs error:       %15.6e\n" $abs_err
  rel_err=$(echo "$abs_err / $de" | bc -l)
  printf "rel error:       %15.6e\n" $rel_err

  echo "==============================================================="
}

if (( $# == 3 ))
then
  # check_stress.sh term sigma_name file.r
  echo "==============================================================="
  echo " check_stress.sh" ${*}
  echo "==============================================================="
  check_stress_term $1 $2 $3
  exit
fi

if (( $# == 1 ))
then
  # check_stress.sh file.r
  echo "==============================================================="
  echo " check_stress.sh" ${*}
  echo "==============================================================="
  check_stress_term "ekin"    "sigma_ekin" $1
  check_stress_term "eps"     "sigma_eps"  $1
  check_stress_term "econf"   "sigma_econf"  $1
  check_stress_term "enl"     "sigma_enl"  $1
  check_stress_term "exc"     "sigma_exc"  $1
  check_stress_term "esr"     "sigma_esr"  $1
  check_stress_term "etotal"  "sigma_eks"  $1
  exit
fi

# incorrect number of arguments. Print usage.
echo " use: $0 [energy_term sigma_name] file.r"
