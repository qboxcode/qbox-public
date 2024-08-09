#!/usr/bin/env python3
# Copyright 2016 The Regents of the University of California
# This file is part of Qbox
#
# qbox_xyz.py: extract sets of atomic positions in xyz format
# from a Qbox output file or from a Qbox sample file using SAX
# incremental parsing
# Optionally fold atoms into the unit cell [-pbc]
#
# use: qbox_xyz.py [-pbc] {file|URL}
import os.path
import xml.sax
import sys
import urllib
from qso import UnitCell
import numpy as np

def usage():
  print ("use: ",sys.argv[0]," [-pbc] {file|URL}")
  sys.exit()

argc=len(sys.argv)
if ( argc < 2 or argc > 3 ):
  usage()

# check if option "-pbc" is used
use_pbc = False
input_source = sys.argv[1]
if ( sys.argv[1] == "-pbc" ):
  if ( argc != 3 ):
    usage()
  use_pbc = True
  input_source = sys.argv[2]

# conversion from Bohr to Angstrom
a0=0.529177

# Qbox output handler to extract and process data
class QboxOutputHandler(xml.sax.handler.ContentHandler):
  def __init__(self):
    self.step = 0
    self.inAtomset = 0
    self.inAtom = 0
    self.inPosition = 0
    self.done_first = False
    self.cell = UnitCell()

  def startElement(self, name, attributes):
    if name == "atomset":
      self.tau=[]
      self.atomname=[]
      self.inAtomset = 1
    elif (name == "unit_cell") & self.inAtomset:
      self.cell.a = [float(s) for s in attributes["a"].split()]
      self.cell.b = [float(s) for s in attributes["b"].split()]
      self.cell.c = [float(s) for s in attributes["c"].split()]
      self.cell.update()
    elif (name == "atom") & self.inAtomset:
      self.atomname.append(attributes["name"])
      self.inAtom = 1
    elif (name == "position") & self.inAtom:
        self.buffer = ""
        self.inPosition = 1

  def characters(self, data):
    if self.inPosition:
      self.buffer += data

  def endElement(self, name):
    if (name == "atom") and self.inAtomset:
      self.inAtom = 0
    if (name == "position") & self.inAtom:
      pos = self.buffer.split()
      x = float(pos[0])
      y = float(pos[1])
      z = float(pos[2])
      r = np.array([x,y,z])
      if ( use_pbc ):
        r=np.array(self.cell.fold_in_ws(r[0],r[1],r[2]))

      self.tau.append([a0*r[0],a0*r[1],a0*r[2]])
      self.inPosition = 0
    elif name == "atomset":
      self.step += 1
      print (len(self.tau))
      av = self.cell.a
      bv = self.cell.b
      cv = self.cell.c
      print (self.step,
      '%.6f'%(a0*(av[0])), '%.6f'%(a0*(av[1])), '%.6f'%(a0*(av[2])),
      '%.6f'%(a0*(bv[0])), '%.6f'%(a0*(bv[1])), '%.6f'%(a0*(bv[2])),
      '%.6f'%(a0*(cv[0])), '%.6f'%(a0*(cv[1])), '%.6f'%(a0*(cv[2])))
      for i in range(len(self.tau)):
        print (self.atomname[i],'%.6f'%self.tau[i][0],\
                                '%.6f'%self.tau[i][1],\
                                '%.6f'%self.tau[i][2])
      self.inAtomset = 0
      self.done_first = True

parser = xml.sax.make_parser()
handler = QboxOutputHandler()
parser.setContentHandler(handler)
# test if input_source is a local file
# if not, process as a URL
if ( os.path.isfile(input_source) ):
  file = open(input_source)
  s = file.read(8192)
  while ( s !="" ):
    parser.feed(s)
    s = file.read(8192)
  file.close()
else:
  # attempt to open as a URL
  try:
    f = urllib.request.urlopen(input_source)
    s = f.read(8192)
    while ( s !="" ):
      parser.feed(s)
      s = f.read(8192)
    f.close()
  except (ValueError,urllib.error.HTTPError) as e:
    print (e)
    sys.exit()

parser.reset()
