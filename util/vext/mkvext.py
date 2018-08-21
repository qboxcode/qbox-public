#!/usr/bin/python
#
# mkvext.py: generate a potential file vext.xml
#
# The function f(x,y,z) defines the potential

import math
import base64
import struct

# example: quadratic potential centered at the origin
def f(x,y,z):
  return 0.5*(x*x+y*y+z*z)

def norm(x):
  return math.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])

a=[16.0,  0.0,  0.0]
b=[ 0.0, 16.0,  0.0]
c=[ 0.0,  0.0, 16.0]

n=[64,64,64]

h=[1.0/n[0],1.0/n[1],1.0/n[2]]

# origin of coordinates in the domain corner
v=[]
for k in range(n[2]):
  if ( k < n[2]/2 ):
    fc = k*h[2]
  else:
    fc = (k-n[2])*h[2]
  for j in range(n[1]):
    if ( j < n[1]/2 ):
      fb = j*h[1]
    else:
      fb = (j-n[1])*h[1]
    for i in range(n[0]):
      if ( i < n[0]/2 ):
        fa = i*h[0]
      else:
        fa = (i-n[0])*h[0]
      x = fa * a[0] + fb * b[0] + fc *c[0]
      y = fa * a[1] + fb * b[1] + fc *c[1]
      z = fa * a[2] + fb * b[2] + fc *c[2]
      v.append(f(float(x),float(y),float(z)))
      #print x,y,z,f(float(x),float(y),float(z))

print '<?xml version="1.0" encoding="UTF-8"?>'
print '<fpmd:function3d xmlns:fpmd="http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0"'
print ' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'
print ' xsi:schemaLocation="http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0 function3d.xsd"'
print ' name="delta_v">'
print '<domain a="'+str(a[0])+" "+str(a[1])+" "+str(a[2])+'"'
print '        b="'+str(b[0])+" "+str(b[1])+" "+str(b[2])+'"'
print '        c="'+str(c[0])+" "+str(c[1])+" "+str(c[2])+'"/>'
print '<grid nx="'+str(n[0])+'" ny="'+str(n[1])+'" nz="'+str(n[2])+'"/>'
print '<grid_function type="double" nx="'+str(n[0])+'" ny="'+str(n[1])+'" nz="'+str(n[2])+'" encoding="base64">'
# encode vector in base64
str="".join(struct.pack('d',t) for t in v)
s=base64.encodestring(str).strip()
print s
print '</grid_function>'
print '</fpmd:function3d>'
