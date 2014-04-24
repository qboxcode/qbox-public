#!/usr/bin/python
# qbox_msd.py: compute mean-square displacement in an MD simulation
# generate plot of <r^2>(t) in gnuplot format
# use: qbox_msd.py species file.r [file.r ...]

import xml.sax
import sys
import math

if len(sys.argv) < 3:
  print "use: ",sys.argv[0]," species file [file ...]"
  sys.exit()

species = sys.argv[1]

# Qbox output handler to extract and process data
class QboxOutputHandler(xml.sax.handler.ContentHandler):
  def __init__(self):
    self.step = 0
    self.readData = 0
    self.readAtom = 0
    self.tau0 = []

  def startElement(self, name, attributes):
    if name == "atomset":
      self.tau=[]
    elif name == "atom":
      self.readAtom = (attributes["species"] == species)
    elif (name == "position") & self.readAtom:
        self.readData = 1
        self.buffer = ""

  def characters(self, data):
    if self.readData:
      self.buffer += data

  def endElement(self, name):
    if (name == "position") & self.readAtom:
      self.readData = 0
      pos = self.buffer.split()
      x = float(pos[0])
      y = float(pos[1])
      z = float(pos[2])
      self.tau.append([x,y,z])
    elif name == "atomset":
      if self.step == 0:
        # copy initial positions to tau0
        self.tau0 = self.tau
      # compute square displacement
      disp2sum = 0.0
      for i in range(len(self.tau)):
        dx = self.tau[i][0]-self.tau0[i][0]
        dy = self.tau[i][1]-self.tau0[i][1]
        dz = self.tau[i][2]-self.tau0[i][2]
        disp2sum += dx*dx + dy*dy + dz*dz
      print '%12.6f'%(disp2sum/len(self.tau))
      self.step += 1

print "# ",species
parser = xml.sax.make_parser()
handler = QboxOutputHandler()
parser.setContentHandler(handler)
for i in range(len(sys.argv)-2):
  infile = sys.argv[i+2]
  parser.parse(infile)

print
print
