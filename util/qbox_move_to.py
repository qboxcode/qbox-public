#!/usr/bin/env python3
# Copyright 2018 The Regents of the University of California
# This file is part of Qbox
#
# qbox_move_to.py: create a qbox input file moving atoms to the positions
# of a given iteration of a simulation, or of a restart file.
#
# use: qbox_move_to.py [iter] {file|URL}
import os.path
import xml.sax
import sys
import urllib
from qso import UnitCell

def usage():
  print ("use: ",sys.argv[0]," [-iter i] {file|URL}")
  sys.exit()

argc=len(sys.argv)
if ( argc != 2 and argc != 4 ):
  usage()

# check if option -iter is used
# "iter" option: extract atomset of iteration iter
# default: iter = 1
iter = 1
iter_option = False
input_source = sys.argv[1]
if ( sys.argv[1] == "-iter" ):
  if ( argc != 4 ):
    usage()
  iter_option = True
  iter = int(sys.argv[2])
  input_source = sys.argv[3]

# Qbox output handler to extract and process data
class QboxOutputHandler(xml.sax.handler.ContentHandler):
  def __init__(self):
    self.step = 0
    self.inAtomset = 0
    self.inAtom = 0
    self.inPosition = 0
    self.inVelocity = 0
    self.done = False
    self.cell = UnitCell()

  def startElement(self, name, attributes):
    if name == "atomset":
      self.tau=[]
      self.vel=[]
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
    elif (name == "velocity") & self.inAtom:
        self.buffer = ""
        self.inVelocity = 1

  def characters(self, data):
    if (self.inPosition) or (self.inVelocity):
      self.buffer += data

  def endElement(self, name):
    if (name == "atom") and self.inAtomset:
      self.inAtom = 0
    if (name == "position") & self.inAtom:
      pos = self.buffer.split()
      self.tau.append([pos[0],pos[1],pos[2]])
      self.inPosition = 0
    if (name == "velocity") & self.inAtom:
      v = self.buffer.split()
      self.vel.append([v[0],v[1],v[2]])
      self.inVelocity = 0
    elif name == "atomset":
      self.step += 1
      if ( self.step == iter ):
        if ( iter_option ):
          print ("#",input_source,"iteration",iter)
        else:
          print ("#",input_source)
        av = self.cell.a
        bv = self.cell.b
        cv = self.cell.c
        print ("set cell ",av[0],av[1],av[2],
               bv[0],bv[1],bv[2],cv[0],cv[1],cv[2])
        for i in range(len(self.tau)):
          print ("move ",self.atomname[i]," to ",\
            self.tau[i][0],self.tau[i][1],self.tau[i][2],\
            self.vel[i][0],self.vel[i][1],self.vel[i][2])
      self.inAtomset = 0
      self.done = ( self.step >= iter )

parser = xml.sax.make_parser()
handler = QboxOutputHandler()
parser.setContentHandler(handler)
# test if input_source is a local file
# if not, process as a URL
if ( os.path.isfile(input_source) ):
  file = open(input_source)
  s = file.read(8192)
  while ( s !="" and not handler.done ):
    parser.feed(s)
    s = file.read(8192)
  file.close()
else:
  # attempt to open as a URL
  try:
    f = urllib.request.urlopen(input_source)
    s = f.read(8192)
    while ( s !="" and not handler.done ):
      parser.feed(s)
      s = f.read(8192)
    f.close()
  except (ValueError,urllib.error.HTTPError) as e:
    print (e)
    sys.exit()

parser.reset()
