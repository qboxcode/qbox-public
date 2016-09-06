#!/usr/bin/python
# Copyright 2016 The Regents of the University of California
# This file is part of Qbox
#
# qbox_maxforce.py: extract largest force from a Qbox output file
# using SAX incremental parsing
#
# use: qbox_maxforce.py {file|URL}
import os.path
import xml.sax
import sys
import urllib2

def usage():
  print "use: ",sys.argv[0]," {file|URL}"
  sys.exit()

argc=len(sys.argv)
if ( argc != 2 ):
  usage()

input_source = sys.argv[1]

# Qbox output handler to extract and process data
class QboxOutputHandler(xml.sax.handler.ContentHandler):
  def __init__(self):
    self.step = 0
    self.inAtomset = 0
    self.inAtom = 0
    self.inForce = 0

  def startElement(self, name, attributes):
    if name == "atomset":
      self.f=[]
      self.atomname=[]
      self.inAtomset = 1
    elif (name == "atom") & self.inAtomset:
      self.atomname.append(attributes["name"])
      self.inAtom = 1
    elif (name == "force") & self.inAtom:
        self.buffer = ""
        self.inForce = 1

  def characters(self, data):
    if self.inForce:
      self.buffer += data

  def endElement(self, name):
    if (name == "atom") and self.inAtomset:
      self.inAtom = 0
    if (name == "force") & self.inAtom:
      force = self.buffer.split()
      fx = float(force[0])
      fy = float(force[1])
      fz = float(force[2])
      self.f.append([fx,fy,fz])
      self.inForce = 0
    elif name == "atomset":
      self.step += 1
      fxmax = 0.0
      fymax = 0.0
      fzmax = 0.0
      x_name = ' '
      y_name = ' '
      z_name = ' '
      for i in range(len(self.f)):
        fx = self.f[i][0]
        if ( fx*fx > fxmax*fxmax ):
          fxmax = fx
          x_name = self.atomname[i]
        fy = self.f[i][1]
        if ( fy*fy > fymax*fymax ):
          fymax = fy
          y_name = self.atomname[i]
        fz = self.f[i][2]
        if ( fz*fz > fzmax*fzmax ):
          fzmax = fz
          z_name = self.atomname[i]
      print '%10.3e'%fxmax,'%-8s'%x_name,\
            '%10.3e'%fymax,'%-8s'%y_name,\
            '%10.3e'%fzmax,'%-8s'%z_name
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
