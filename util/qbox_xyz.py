#!/usr/bin/python
# Copyright 2016 The Regents of the University of California
# This file is part of Qbox
#
# qbox_xyz.py: extract sets of atomic positions in xyz format
# from a Qbox output file or from a Qbox sample file using SAX
# incremental parsing
#
# use: qbox_xyz.py [-first] {file|URL}
import os.path
import xml.sax
import sys
import urllib2

def usage():
  print "use: ",sys.argv[0]," [-first] {file|URL}"
  sys.exit()

argc=len(sys.argv)
if ( argc < 2 or argc > 3 ):
  usage()

# check if option "-first" is used
# "-first" option: extract first atomset only
# default: extract all atomsets
first_only = False
input_source = sys.argv[1]
if ( sys.argv[1] == "-first" ):
  if ( argc != 3 ):
    usage()
  first_only = True
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
      x = a0*float(pos[0])
      y = a0*float(pos[1])
      z = a0*float(pos[2])
      self.tau.append([x,y,z])
      self.inPosition = 0
    elif name == "atomset":
      self.step += 1
      print len(self.tau)
      avec = self.cell_a.split()
      bvec = self.cell_b.split()
      cvec = self.cell_c.split()
      print self.step,\
      '%.6f'%(a0*float(avec[0])),\
      '%.6f'%(a0*float(avec[1])),\
      '%.6f'%(a0*float(avec[2])),\
      '%.6f'%(a0*float(bvec[0])),\
      '%.6f'%(a0*float(bvec[1])),\
      '%.6f'%(a0*float(bvec[2])),\
      '%.6f'%(a0*float(cvec[0])),\
      '%.6f'%(a0*float(cvec[1])),\
      '%.6f'%(a0*float(cvec[2]))
      for i in range(len(self.tau)):
        print self.atomname[i],'%.6f'%self.tau[i][0],\
                               '%.6f'%self.tau[i][1],\
                               '%.6f'%self.tau[i][2]
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
  while ( s !="" and not (first_only and handler.done_first) ):
    parser.feed(s)
    s = file.read(8192)
  file.close()
else:
  # attempt to open as a URL
  try:
    f = urllib2.urlopen(input_source)
    s = f.read(8192)
    while ( s !="" and not (first_only and handler.done_first) ):
      parser.feed(s)
      s = f.read(8192)
    f.close()
  except (ValueError,urllib2.HTTPError) as e:
    print e
    sys.exit()

parser.reset()
