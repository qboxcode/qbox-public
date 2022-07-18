#!/usr/bin/env python2
# qbox_angle.py
# extract angle defined by three atoms from Qbox output
# use: qbox_angle.py name1 name2 name3 file.r
import xml.sax
import sys
import math

if len(sys.argv) != 5:
  print "use: ",sys.argv[0]," name1 name2 name3 file.r"
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

  def startElement(self, name, attributes):
    if name == "atomset":
      self.buffer1 = ""
      self.buffer2 = ""
      self.buffer3 = ""
      self.atom1done = 0
      self.atom2done = 0
      self.atom3done = 0
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
      r1 = (float(pos1[0]),float(pos1[1]),float(pos1[2]))
      pos2 = self.buffer2.split()
      r2 = (float(pos2[0]),float(pos2[1]),float(pos2[2]))
      pos3 = self.buffer3.split()
      r3 = (float(pos3[0]),float(pos3[1]),float(pos3[2]))
      #print "r1: ",r1
      #print "r2: ",r2
      #print "r3: ",r3
      # e12 = normalized r12
      r12 = (r1[0]-r2[0],r1[1]-r2[1],r1[2]-r2[2])
      fac12 = 1.0/norm(r12)
      e12 = (fac12*r12[0],fac12*r12[1],fac12*r12[2])
      # e32 = normalized r32
      r32 = (r3[0]-r2[0],r3[1]-r2[1],r3[2]-r2[2])
      fac32 = 1.0/norm(r32)
      e32 = (fac32*r32[0],fac32*r32[1],fac32*r32[2])
      sp = e12[0]*e32[0] + e12[1]*e32[1] + e12[2]*e32[2]
      c = sp
      c = max(-1.0,min(1.0,sp))
      a = (180.0/math.pi)*math.acos(c)
      print '%.4f' % a
    elif name == "position":
      self.readPos1 = 0
      self.readPos2 = 0
      self.readPos3 = 0

name1 = sys.argv[1]
name2 = sys.argv[2]
name3 = sys.argv[3]
filename = sys.argv[4]
parser = xml.sax.make_parser()
handler = QboxOutputHandler()
parser.setContentHandler(handler)
parser.parse(sys.argv[4])
