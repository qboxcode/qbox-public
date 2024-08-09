#!/usr/bin/env python3
# qbox_torsion.py
# extract torsion angle defined by four atoms from Qbox output
# use: qbox_torsion.py name1 name2 name3 name4 file.r
import xml.sax
import sys
import math
from qso import UnitCell
import numpy as np

use_msg = "use: "+sys.argv[0]+" [-pbc] name1 name2 name3 name4 file.r"
if len(sys.argv) < 6 or len(sys.argv) > 7:
  print(use_msg)
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
    self.inAtomSet = False
    self.cell = UnitCell()

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
      self.inAtomSet = True
    elif (name == "unit_cell") and self.inAtomSet:
      self.cell.a = [float(s) for s in attributes["a"].split()]
      self.cell.b = [float(s) for s in attributes["b"].split()]
      self.cell.c = [float(s) for s in attributes["c"].split()]
      self.cell.update()
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
      r1 = np.array([float(pos1[0]),float(pos1[1]),float(pos1[2])])
      pos2 = self.buffer2.split()
      r2 = np.array([float(pos2[0]),float(pos2[1]),float(pos2[2])])
      pos3 = self.buffer3.split()
      r3 = np.array([float(pos3[0]),float(pos3[1]),float(pos3[2])])
      pos4 = self.buffer4.split()
      r4 = np.array([float(pos4[0]),float(pos4[1]),float(pos4[2])])

      # e12 = normalized r12
      r12 = r1 - r2
      if ( use_pbc ):
        r12=np.array(self.cell.fold_in_ws(r12[0],r12[1],r12[2]))
      e12 = r12 / norm(r12)
      # e32 = normalized r32
      r32 = r3 - r2
      if ( use_pbc ):
        r32=np.array(self.cell.fold_in_ws(r32[0],r32[1],r32[2]))
      e32 = r32 / norm(r32)
      # e43 = normalized r43
      r43 = r4 - r3
      if ( use_pbc ):
        r43=np.array(self.cell.fold_in_ws(r43[0],r43[1],r43[2]))
      e43 = r43 / norm(r43)
      # e23 = - e32
      e23 = -e32

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

      print ('%.4f' % a)
      self.inAtomSet = True
    elif name == "position":
      self.readPos1 = 0
      self.readPos2 = 0
      self.readPos3 = 0
      self.readPos4 = 0

if ( len(sys.argv) == 7 ):
  if ( sys.argv[1] == "-pbc" ):
    use_pbc = True
    name1 = sys.argv[2]
    name2 = sys.argv[3]
    name3 = sys.argv[4]
    name4 = sys.argv[5]
    filename = sys.argv[6]
  else:
    print(use_msg)
    sys.exit()
else:
  use_pbc = False
  name1 = sys.argv[1]
  name2 = sys.argv[2]
  name3 = sys.argv[3]
  name4 = sys.argv[4]
  filename = sys.argv[5]

parser = xml.sax.make_parser()
handler = QboxOutputHandler()
parser.setContentHandler(handler)
parser.parse(filename)
