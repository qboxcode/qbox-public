#!/usr/bin/env python
# generate a Qbox input file from an xyz file
# use: xyz2qbox.py xyzfile
import sys
f = open(sys.argv[1])
line = f.readline()
buf = line.split()
nat = int(buf[0])
line = f.readline()
print "#",line,
for i in range(nat):
  line = f.readline()
  buf = line.split()
  name=buf[0]
  x = float(buf[1])/0.529177
  y = float(buf[2])/0.529177
  z = float(buf[3])/0.529177
  print "atom ",name+str(i+1)," ",name+"_species",'%9.4f'%x,'%9.4f'%y,'%9.4f'%z
