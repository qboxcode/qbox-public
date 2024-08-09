#!/usr/bin/env python3
# qbox_distance.py
# extract distance between two atoms from Qbox output
# use: qbox_distance.py name1 name2 file.r
import xml.sax
import sys
import math
from qso import UnitCell
import numpy as np

use_msg = "use: "+sys.argv[0]+" [-pbc] name1 name2 file.r"
if len(sys.argv) < 4 or len(sys.argv) > 5:
  print(use_msg)
  sys.exit()

name1 = ""
name2 = ""

# Qbox output handler to extract and process <atomset>
class QboxOutputHandler(xml.sax.handler.ContentHandler):
  def __init__(self):
    self.readPos1 = 0
    self.readPos2 = 0
    self.buffer1 = ""
    self.buffer2 = ""
    self.inAtomSet = False
    self.cell = UnitCell()

  def startElement(self, name, attributes):
    if name == "atomset":
      self.buffer1 = ""
      self.buffer2 = ""
      self.atom1done = 0
      self.atom2done = 0
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
    elif name == "position":
      self.buffer = ""

  def characters(self, data):
    if self.readPos1:
      self.buffer1 += data
    elif self.readPos2:
      self.buffer2 += data

  def endElement(self, name):
    if name == "atomset":
      pos1 = self.buffer1.split()
      pos2 = self.buffer2.split()
      r1x = float(pos1[0])
      r1y = float(pos1[1])
      r1z = float(pos1[2])
      r2x = float(pos2[0])
      r2y = float(pos2[1])
      r2z = float(pos2[2])
      r1 = np.array([r1x,r1y,r1z])
      r2 = np.array([r2x,r2y,r2z])
      dr = r2 - r1
      if ( use_pbc ):
        dr=np.array(self.cell.fold_in_ws(dr[0],dr[1],dr[2]))
      print ('%.4f' % math.sqrt(dr[0]**2+dr[1]**2+dr[2]**2))
      self.inAtomSet = False
    elif name == "position":
      self.readPos1 = 0
      self.readPos2 = 0

if ( len(sys.argv) == 5 ):
  if ( sys.argv[1] == "-pbc" ):
    use_pbc = True
    name1 = sys.argv[2]
    name2 = sys.argv[3]
    filename = sys.argv[4]
  else:
    print(use_msg)
    sys.exit()
else:
  use_pbc = False
  name1 = sys.argv[1]
  name2 = sys.argv[2]
  filename = sys.argv[3]

parser = xml.sax.make_parser()
handler = QboxOutputHandler()
parser.setContentHandler(handler)
parser.parse(filename)
