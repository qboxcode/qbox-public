#!/usr/bin/python
# qbox_repair_h2o.py: repair broken h2o molecules in a Qbox sys file
# move hydrogen atoms across periodic boundaries to repair molecules
# use: qbox_repair_h2o.py file.sys

import sys
import math

olist = []
hlist = []

def distance(a,b,sx,sy,sz):
  return math.sqrt((a[3]-b[3]-sx)**2+(a[4]-b[4]-sy)**2+(a[5]-b[5]-sz)**2)

def fold_in_ws(atom):
  x = atom[3]
  y = atom[4]
  z = atom[5]
  while x > 0.5*a_cell + 1.e-5:
    x -= a_cell
  while x < -0.5*a_cell - 1.e-5:
    x += a_cell
  while y > 0.5*b_cell + 1.e-5:
    y -= b_cell
  while y < -0.5*b_cell - 1.e-5:
    y += b_cell
  while z > 0.5*c_cell + 1.e-5:
    z -= c_cell
  while z < -0.5*c_cell - 1.e-5:
    z += c_cell
  atom[3] = x
  atom[4] = y
  atom[5] = z

f = open(sys.argv[1])
for line in f:
  l = line.split()
  if ( l[0] == "set" ) & ( l[1] == "cell" ):
    a_cell = float(l[2])
    b_cell = float(l[6])
    c_cell = float(l[10])
    print line
  elif ( l[0] == "species" ):
    print line
  elif ( l[0] == "atom" ) & ( l[2] == "oxygen" ):
    olist.append([l[0],l[1],l[2],float(l[3]),float(l[4]),float(l[5])])
  elif ( l[0] == "atom" ) & ( l[2] == "hydrogen" ):
    hlist.append([l[0],l[1],l[2],float(l[3]),float(l[4]),float(l[5])])

for o in olist:
  fold_in_ws(o)
for h in hlist:
  fold_in_ws(h)

for h in hlist:
  # find nearest oxygen atom in olist
  mindist = 1.e10
  for o in olist:
    # compute minimal distance o-h
    for sx in [-a_cell,0,a_cell]:
      for sy in [-b_cell,0,b_cell]:
        for sz in [-c_cell,0,c_cell]:
          d = distance(o,h,sx,sy,sz)
          # print "dist(",sx,sy,sz,") = ",d
          if d < mindist:
            mindist = d
            sx_min = sx
            sy_min = sy
            sz_min = sz
            nearest_o = o;
  # print "min shift is: ",sx_min,sy_min,sz_min
  if ( sx_min != 0 ) | ( sy_min != 0 ) | ( sz_min != 0 ):
    print "# current ",h[1]," at ", h[3],h[4],h[5]
    print "# nearest O is at ", nearest_o[3],nearest_o[4],nearest_o[5]
    print "# move ",h[1]," by ", sx_min, sy_min, sz_min
    h[3] += sx_min
    h[4] += sy_min
    h[5] += sz_min

for o in olist:
  print o[0],o[1],o[2],'%10.5f'%o[3],'%10.5f'%o[4],'%10.5f'%o[5]
for h in hlist:
  print h[0],h[1],h[2],'%10.5f'%h[3],'%10.5f'%h[4],'%10.5f'%h[5]
