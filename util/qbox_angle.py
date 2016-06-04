#!/usr/bin/python
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
      self.pos1 = self.buffer1.split()
      self.pos2 = self.buffer2.split()
      self.pos3 = self.buffer3.split()
      self.r1x = float(self.pos1[0])
      self.r1y = float(self.pos1[1])
      self.r1z = float(self.pos1[2])
      self.r2x = float(self.pos2[0])
      self.r2y = float(self.pos2[1])
      self.r2z = float(self.pos2[2])
      self.r3x = float(self.pos3[0])
      self.r3y = float(self.pos3[1])
      self.r3z = float(self.pos3[2])
      #print "r1: ",self.r1x,self.r1y,self.r1z
      #print "r2: ",self.r2x,self.r2y,self.r2z
      #print "r3: ",self.r3x,self.r3y,self.r3z
      # normalized r12
      r12x = self.r1x - self.r2x
      r12y = self.r1y - self.r2y
      r12z = self.r1z - self.r2z
      fac12 = 1.0/math.sqrt(r12x*r12x+r12y*r12y+r12z*r12z)
      r12x = fac12 * r12x
      r12y = fac12 * r12y
      r12z = fac12 * r12z
      # normalized r32
      r32x = self.r3x - self.r2x
      r32y = self.r3y - self.r2y
      r32z = self.r3z - self.r2z
      fac32 = 1.0/math.sqrt(r32x*r32x+r32y*r32y+r32z*r32z)
      r32x = fac32 * r32x
      r32y = fac32 * r32y
      r32z = fac32 * r32z
      sp = r12x*r32x + r12y*r32y + r12z*r32z
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
