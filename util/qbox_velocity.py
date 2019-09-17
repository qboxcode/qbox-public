#!/usr/bin/python
# qbox_velocity.py
# extract velocity of an atom from Qbox output
# use: qbox_velocity.py atom_name file.r
import xml.sax
import sys
import math

if len(sys.argv) != 3:
  print "use: ",sys.argv[0]," atom_name file.r"
  sys.exit()

# Qbox output handler to extract and process <atomset>
class QboxOutputHandler(xml.sax.handler.ContentHandler):
  def __init__(self):
    self.readVelocity = 0
    self.buffer = ""

  def startElement(self, name, attributes):
    if name == "atomset":
      self.buffer = ""
    elif name == "atom":
      if attributes["name"] == atom_name:
        self.readVelocity = 1
    elif name == "velocity":
      self.buffer = ""

  def characters(self, data):
    if self.readVelocity:
      self.buffer += data

  def endElement(self, name):
    if name == "velocity":
      if self.readVelocity:
        velocity = self.buffer.split()
        vx = float(velocity[0])
        vy = float(velocity[1])
        vz = float(velocity[2])
        print '%.8f'%vx,'%.8f'%vy,'%.8f'%vz
        self.readVelocity = 0

atom_name = sys.argv[1]
filename = sys.argv[2]
parser = xml.sax.make_parser()
handler = QboxOutputHandler()
parser.setContentHandler(handler)
parser.parse(sys.argv[2])
