#!/usr/bin/env python3
# qbox_angle.py
# extract angle defined by three atoms from Qbox output
# use: qbox_angle.py name1 name2 name3 file.r
import xml.sax
import sys
import math
from qso import UnitCell
import numpy as np

use_msg = "use: "+sys.argv[0]+" [-pbc] name1 name2 name3 file.r"
if len(sys.argv) < 5 or len(sys.argv) > 6:
  print(use_msg)
  sys.exit()

name1 = ""
name2 = ""
name3 = ""

def norm(x):
  return math.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])

# Qbox output handler to extract and process <atomset>
class QboxOutputHandler(xml.sax.handler.ContentHandler):
  def __init__(self):
    self.readPos1 = 0
    self.readPos2 = 0
    self.readPos3 = 0
    self.buffer1 = ""
    self.buffer2 = ""
    self.buffer3 = ""
    self.inAtomSet = False
    self.cell = UnitCell()

  def startElement(self, name, attributes):
    if name == "atomset":
      self.buffer1 = ""
      self.buffer2 = ""
      self.buffer3 = ""
      self.atom1done = 0
      self.atom2done = 0
      self.atom3done = 0
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
    elif name == "position":
      self.buffer = ""

  def characters(self, data):
    if self.readPos1:
      self.buffer1 += data
    elif self.readPos2:
      self.buffer2 += data
    elif self.readPos3:
      self.buffer3 += data

  def endElement(self, name):
    if name == "atomset":
      pos1 = self.buffer1.split()
      r1 = np.array([float(pos1[0]),float(pos1[1]),float(pos1[2])])
      pos2 = self.buffer2.split()
      r2 = np.array([float(pos2[0]),float(pos2[1]),float(pos2[2])])
      pos3 = self.buffer3.split()
      r3 = np.array([float(pos3[0]),float(pos3[1]),float(pos3[2])])
      # e12 = normalized r12
      r12 = r1 - r2
      if ( use_pbc ):
        r12=np.array(self.cell.fold_in_ws(r12[0],r12[1],r12[2]))
      e12 = r12 / norm(r12)
      print("r12=",r12)
      # e32 = normalized r32
      r32 = r3 - r2
      if ( use_pbc ):
        r32=np.array(self.cell.fold_in_ws(r32[0],r32[1],r32[2]))
      e32 = r32 / norm(r32)
      print("r32=",r32)
      sp = np.dot(e12,e32)
      c = sp
      c = max(-1.0,min(1.0,sp))
      a = (180.0/math.pi)*math.acos(c)
      print ('%.4f' % a)
      self.inAtomSet = True
    elif name == "position":
      self.readPos1 = 0
      self.readPos2 = 0
      self.readPos3 = 0

if ( len(sys.argv) == 6 ):
  if ( sys.argv[1] == "-pbc" ):
    use_pbc = True
    name1 = sys.argv[2]
    name2 = sys.argv[3]
    name3 = sys.argv[4]
    filename = sys.argv[5]
  else:
    print(use_msg)
    sys.exit()
else:
  use_pbc = False
  name1 = sys.argv[1]
  name2 = sys.argv[2]
  name3 = sys.argv[3]
  filename = sys.argv[4]

parser = xml.sax.make_parser()
handler = QboxOutputHandler()
parser.setContentHandler(handler)
parser.parse(filename)
