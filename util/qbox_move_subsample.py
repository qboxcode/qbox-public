#!/usr/bin/python
# Copyright 2018 The Regents of the University of California
# This file is part of Qbox
#
# qbox_move_subsample.py: create a qbox input file moving atoms to the positions
# of iterations of a simulation subsampled every "interval" number of steps
# Each set of moves is followed by a command cmd given as an argument
#
# use: qbox_move_subsample.py interval cmd {file|URL}
import os.path
import xml.sax
import sys
import urllib2

def usage():
  print "use: ",sys.argv[0]," interval cmd {file|URL}"
  sys.exit()

argc=len(sys.argv)
if ( argc != 4 ):
  usage()

interval = int(sys.argv[1])
cmd = sys.argv[2]
input_source = sys.argv[3]

# Qbox output handler to extract and process data
class QboxOutputHandler(xml.sax.handler.ContentHandler):
  def __init__(self):
    self.step = 0
    self.inAtomset = 0
    self.inAtom = 0
    self.inPosition = 0

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
      if ( self.step % interval == 0 ):
        print "#",input_source,"iteration",self.step
        avec = self.cell_a.split()
        bvec = self.cell_b.split()
        cvec = self.cell_c.split()
        print "set cell ",avec[0],avec[1],avec[2],\
          bvec[0],bvec[1],bvec[2],\
          cvec[0],cvec[1],cvec[2]
        for i in range(len(self.tau)):
          print "move ",self.atomname[i]," to ",\
            self.tau[i][0],self.tau[i][1],self.tau[i][2]
        print cmd
      self.inAtomset = 0

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
    f = urllib2.urlopen(input_source)
    s = f.read(8192)
    while ( s !="" ):
      parser.feed(s)
      s = f.read(8192)
    f.close()
  except (ValueError,urllib2.HTTPError) as e:
    print e
    sys.exit()

parser.reset()
