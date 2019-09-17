#!/usr/bin/python
# qbox_force.py
# extract force of an atom from Qbox output
# use: qbox_force.py atom_name file.r
import xml.sax
import sys
import math

if len(sys.argv) != 3:
  print "use: ",sys.argv[0]," atom_name file.r"
  sys.exit()

# Qbox output handler to extract and process <atomset>
class QboxOutputHandler(xml.sax.handler.ContentHandler):
  def __init__(self):
    self.readForce = 0
    self.buffer = ""

  def startElement(self, name, attributes):
    if name == "atomset":
      self.buffer = ""
    elif name == "atom":
      if attributes["name"] == atom_name:
        self.readForce = 1
    elif name == "force":
      self.buffer = ""

  def characters(self, data):
    if self.readForce:
      self.buffer += data

  def endElement(self, name):
    if name == "force":
      if self.readForce:
        force = self.buffer.split()
        fx = float(force[0])
        fy = float(force[1])
        fz = float(force[2])
        print '%.8f'%fx,'%.8f'%fy,'%.8f'%fz
        self.readForce = 0

atom_name = sys.argv[1]
filename = sys.argv[2]
parser = xml.sax.make_parser()
handler = QboxOutputHandler()
parser.setContentHandler(handler)
parser.parse(sys.argv[2])
