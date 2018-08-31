#!/usr/bin/python
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
import urllib2

def usage():
  print "use: ",sys.argv[0]," [-iter i] {file|URL}"
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
    self.done = False

  def startElement(self, name, attributes):
    if name == "atomset":
      self.tau=[]
      self.atomname=[]
      self.inAtomset = 1
    elif (name == "unit_cell") & self.inAtomset:
      self.cell_a = attributes["a"]
      self.cell_b = attributes["b"]
      self.cell_c = attributes["c"]
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
      self.tau.append([pos[0],pos[1],pos[2]])
      self.inPosition = 0
    elif name == "atomset":
      self.step += 1
      if ( self.step == iter ):
        if ( iter_option ):
          print "#",input_source,"iteration",iter
        else:
          print "#",input_source
        avec = self.cell_a.split()
        bvec = self.cell_b.split()
        cvec = self.cell_c.split()
        print "set cell ",avec[0],avec[1],avec[2],\
          bvec[0],bvec[1],bvec[2],\
          cvec[0],cvec[1],cvec[2]
        for i in range(len(self.tau)):
          print "move ",self.atomname[i]," to ",\
            self.tau[i][0],self.tau[i][1],self.tau[i][2]
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
    f = urllib2.urlopen(input_source)
    s = f.read(8192)
    while ( s !="" and not handler.done ):
      parser.feed(s)
      s = f.read(8192)
    f.close()
  except (ValueError,urllib2.HTTPError) as e:
    print e
    sys.exit()

parser.reset()
