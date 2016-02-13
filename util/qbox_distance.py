#!/usr/bin/python
# qbox_distance.py
# extract distance between two atoms from Qbox output
# use: qbox_distance.py name1 name2 file.r
import xml.sax
import sys
import math

if len(sys.argv) != 4:
  print "use: ",sys.argv[0]," name1 name2 file.r"
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

  def startElement(self, name, attributes):
    if name == "atomset":
      self.buffer1 = ""
      self.buffer2 = ""
      self.atom1done = 0
      self.atom2done = 0
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
      self.pos1 = self.buffer1.split()
      self.pos2 = self.buffer2.split()
      self.r1x = float(self.pos1[0])
      self.r1y = float(self.pos1[1])
      self.r1z = float(self.pos1[2])
      self.r2x = float(self.pos2[0])
      self.r2y = float(self.pos2[1])
      self.r2z = float(self.pos2[2])
      #print "r1: ",self.r1x,self.r1y,self.r1z
      #print "r2: ",self.r2x,self.r2y,self.r2z
      print '%.4f' % math.sqrt((self.r1x-self.r2x)**2+
                               (self.r1y-self.r2y)**2+
                               (self.r1z-self.r2z)**2)
    elif name == "position":
      self.readPos1 = 0
      self.readPos2 = 0

name1 = sys.argv[1]
name2 = sys.argv[2]
filename = sys.argv[3]
parser = xml.sax.make_parser()
handler = QboxOutputHandler()
parser.setContentHandler(handler)
parser.parse(sys.argv[3])
