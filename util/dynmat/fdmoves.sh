#!/bin/bash
# use: fdmoves.sh atom1 atom2  .. atomN
# Note: atom names must appear in the same order as in the Qbox <atomset>
#
# generate moves for finite difference calculations of the dynamical matrix

# h: atomic displacement in (a.u.)
h=0.01
# nitscf: number of scf iterations for each displacement
nitscf=20

for atom in ${*}
do
  # moves along x direction
  echo move $atom by  $h 0 0
  echo run 0 $nitscf
  echo move $atom by -$h 0 0

  echo move $atom by -$h 0 0
  echo run 0 $nitscf
  echo move $atom by  $h 0 0

  # moves along y direction
  echo move $atom by  0  $h  0
  echo run 0 $nitscf
  echo move $atom by  0 -$h  0

  echo move $atom by  0 -$h  0
  echo run 0 $nitscf
  echo move $atom by  0  $h  0

  # moves along z direction
  echo move $atom by  0  0  $h
  echo run 0 $nitscf
  echo move $atom by  0  0 -$h

  echo move $atom by  0  0 -$h
  echo run 0 $nitscf
  echo move $atom by  0  0  $h

done
