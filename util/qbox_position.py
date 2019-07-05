#!/usr/bin/python
# qbox_position.py
# extract position of an atom from Qbox output
# use: qbox_position.py atom_name file.r
import xml.sax
import sys
import math

if len(sys.argv) != 3:
  print "use: ",sys.argv[0]," atom_name file.r"
  sys.exit()

# Qbox output handler to extract and process <atomset>
class QboxOutputHandler(xml.sax.handler.ContentHandler):
  def __init__(self):
    self.readPos = 0
    self.buffer = ""

  def startElement(self, name, attributes):
    if name == "atomset":
      self.buffer = ""
    elif name == "atom":
      if attributes["name"] == atom_name:
        self.readPos = 1
    elif name == "position":
      self.buffer = ""

  def characters(self, data):
    if self.readPos:
      self.buffer += data

  def endElement(self, name):
    if name == "position":
      if self.readPos:
        pos = self.buffer.split()
        rx = float(pos[0])
        ry = float(pos[1])
        rz = float(pos[2])
        print '%.8f'%rx,'%.8f'%ry,'%.8f'%rz
        self.readPos = 0

atom_name = sys.argv[1]
filename = sys.argv[2]
parser = xml.sax.make_parser()
handler = QboxOutputHandler()
parser.setContentHandler(handler)
parser.parse(sys.argv[2])
