#!/bin/bash
exe=../../src/qb
mpirun -np 1 $exe -server qbin_0 qbout_0 > qblog_0 &
mpirun -np 1 $exe -server qbin_1 qbout_1 > qblog_1 &
./twin qbin qbout > twin.out &
wait
