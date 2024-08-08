#!/usr/bin/env python3
# fix discontinuities in dipole data
# use: fix_dipole.py ax ay az bx by bz cx cy cz dipole.dat
# dipole.dat contains the sequence of dipole values, 3 numbers per line
import sys
import numpy as np
from qso import UnitCell

assert(len(sys.argv)==11)
ax=float(sys.argv[1])
ay=float(sys.argv[2])
az=float(sys.argv[3])
bx=float(sys.argv[4])
by=float(sys.argv[5])
bz=float(sys.argv[6])
cx=float(sys.argv[7])
cy=float(sys.argv[8])
cz=float(sys.argv[9])
cell=UnitCell(ax,ay,az,bx,by,bz,cx,cy,cz)

f = open(sys.argv[10],'r')
n = 0
for line in f:
  n += 1
  buf = line.split()
  if ( n == 1 ):
    xm = float(buf[0])
    ym = float(buf[1])
    zm = float(buf[2])
    rm = np.array([xm,ym,zm])
    print ('%14.8f'%rm[0],'%14.8f'%rm[1],'%14.8f'%rm[2])
  else:
    x = float(buf[0])
    y = float(buf[1])
    z = float(buf[2])
    r = np.array([x,y,z])
    dr = r - rm
    dr=np.array(cell.fold_in_ws(dr[0],dr[1],dr[2]))
    r = rm + dr
    print ('%14.8f'%r[0],'%14.8f'%r[1],'%14.8f'%r[2])
    rm = r
