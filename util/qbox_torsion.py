#!/usr/bin/python
# qbox_torsion.py
# extract torsion angle defined by four atoms from Qbox output
# use: qbox_torsion.py name1 name2 name3 name4 file.r
import xml.sax
import sys
import math

if len(sys.argv) != 6:
  print "use: ",sys.argv[0]," name1 name2 name3 name4 file.r"
  sys.exit()

name1 = ""
name2 = ""
name3 = ""
name4 = ""

def norm(x):
  return math.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])

def cross_product(a,b):
# cross product c = a ^ b
  cx = a[1] * b[2] - a[2] * b[1]
  cy = a[2] * b[0] - a[0] * b[2]
  cz = a[0] * b[1] - a[1] * b[0]
  return (cx,cy,cz)

# Qbox output handler to extract and process <atomset>
class QboxOutputHandler(xml.sax.handler.ContentHandler):
  def __init__(self):
    self.readPos1 = 0
    self.readPos2 = 0
    self.readPos3 = 0
    self.readPos4 = 0
    self.buffer1 = ""
    self.buffer2 = ""
    self.buffer3 = ""
    self.buffer4 = ""

  def startElement(self, name, attributes):
    if name == "atomset":
      self.buffer1 = ""
      self.buffer2 = ""
      self.buffer3 = ""
      self.buffer4 = ""
      self.atom1done = 0
      self.atom2done = 0
      self.atom3done = 0
      self.atom4done = 0
    elif name == "atom":
      self.atom_name = attributes["name"]
      if self.atom_name == name1:
        self.readPos1 = 1
      elif self.atom_name == name2:
        self.readPos2 = 1
      elif self.atom_name == name3:
        self.readPos3 = 1
      elif self.atom_name == name4:
        self.readPos4 = 1
    elif name == "position":
      self.buffer = ""

  def characters(self, data):
    if self.readPos1:
      self.buffer1 += data
    elif self.readPos2:
      self.buffer2 += data
    elif self.readPos3:
      self.buffer3 += data
    elif self.readPos4:
      self.buffer4 += data

  def endElement(self, name):
    if name == "atomset":
      pos1 = self.buffer1.split()
      r1 = (float(pos1[0]),float(pos1[1]),float(pos1[2]))
      pos2 = self.buffer2.split()
      r2 = (float(pos2[0]),float(pos2[1]),float(pos2[2]))
      pos3 = self.buffer3.split()
      r3 = (float(pos3[0]),float(pos3[1]),float(pos3[2]))
      pos4 = self.buffer4.split()
      r4 = (float(pos4[0]),float(pos4[1]),float(pos4[2]))

      #print "r1: ",r1
      #print "r2: ",r2
      #print "r3: ",r3
      #print "r4: ",r4
      # e12 = normalized r12
      r12 = (r1[0]-r2[0],r1[1]-r2[1],r1[2]-r2[2])
      fac12 = 1.0/norm(r12)
      e12 = (fac12*r12[0],fac12*r12[1],fac12*r12[2])
      # e32 = normalized r32
      r32 = (r3[0]-r2[0],r3[1]-r2[1],r3[2]-r2[2])
      fac32 = 1.0/norm(r32)
      e32 = (fac32*r32[0],fac32*r32[1],fac32*r32[2])
      # e43 = normalized r43
      r43 = (r4[0]-r3[0],r4[1]-r3[1],r4[2]-r3[2])
      fac43 = 1.0/norm(r43)
      e43 = (fac43*r43[0],fac43*r43[1],fac43*r43[2])
      # e23 = - e32
      e23 = (-e32[0],-e32[1],-e32[2])

      u = cross_product(e12,e32)
      v = cross_product(e23,e43)

      a = 0.0
      unorm = norm(u)
      vnorm = norm(v)
      if (unorm!=0.0) and (vnorm!=0.0):
        u = (u[0]/unorm,u[1]/unorm,u[2]/unorm)
        v = (v[0]/vnorm,v[1]/vnorm,v[2]/vnorm)
        uv = u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
        cc = max(min(uv,1.0),-1.0)

        w = cross_product(u,v)
        we32 = e32[0]*w[0] + e32[1]*w[1] + e32[2]*w[2]
        ss = max(min(we32,1.0),-1.0)
        a = (180.0/math.pi) * math.atan2(ss,cc)

      print '%.4f' % a
    elif name == "position":
      self.readPos1 = 0
      self.readPos2 = 0
      self.readPos3 = 0
      self.readPos4 = 0

name1 = sys.argv[1]
name2 = sys.argv[2]
name3 = sys.argv[3]
name4 = sys.argv[4]
filename = sys.argv[5]
parser = xml.sax.make_parser()
handler = QboxOutputHandler()
parser.setContentHandler(handler)
parser.parse(sys.argv[5])
