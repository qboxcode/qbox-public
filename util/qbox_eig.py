#!/usr/bin/python
# qbox_eig.py: extract eigenvalues from Qbox output
# use: qbox_eig.py kpoint n ispin file.r
# extracts eigenvalue n at (ispin,kpoint)
# note: ispin = 0..1, first eigenvalue is n=1

import xml.sax
import sys
import math

argc=len(sys.argv)
if ( not ( argc in [3,4,6,7] ) ):
  print "use: ",sys.argv[0]," [ispin] [kx ky kz] n file.r"
  print " ispin = 0..1, n = 1..neig"
  sys.exit()

if argc == 7:
  # ispin kx ky kz n file.r
  ispin = int(sys.argv[1])
  kx = float(sys.argv[2])
  ky = float(sys.argv[3])
  kz = float(sys.argv[4])
  n = int(sys.argv[5])
  infile = sys.argv[6]
elif argc == 6:
  # kx ky kz n file.r
  ispin = 0
  kx = float(sys.argv[1])
  ky = float(sys.argv[2])
  kz = float(sys.argv[3])
  n = int(sys.argv[4])
  infile = sys.argv[5]
elif argc == 4:
  # ispin n file.r
  ispin = int(sys.argv[1])
  kx = 0.0
  ky = 0.0
  kz = 0.0
  n = int(sys.argv[2])
  infile = sys.argv[3]
elif argc == 3:
  ispin = 0
  kx = 0.0
  ky = 0.0
  kz = 0.0
  n = int(sys.argv[1])
  infile = sys.argv[2]

print "# ",infile," ispin=",ispin, " n=", n, " k=", kx, ky, kz

# Qbox output handler to extract and process data
class QboxOutputHandler(xml.sax.handler.ContentHandler):
  def __init__(self):
    self.iter = 0
    self.readData = 0

  def startElement(self, name, attributes):
    if name == "eigenvalues":
      self.n = attributes["n"]
      self.spin = attributes["spin"]
      self.kpoint = attributes["kpoint"]
      self.readData = 1
      self.iter += 1
      self.buffer = ""

  def characters(self, data):
    if self.readData:
      self.buffer += data

  def endElement(self, name):
    if name == "eigenvalues":
      self.readData = 0
      isp = int(self.spin)
      self.kp = self.kpoint.split()
      dx = kx-float(self.kp[0])
      dy = ky-float(self.kp[1])
      dz = kz-float(self.kp[2])
      if isp == ispin and math.sqrt(dx*dx+dy*dy+dz*dz) < 1.e-6:
        self.print_eig()

  def print_eig(self):
    self.e = self.buffer.split()
    if n > int(self.n):
      print "n>neig: neig=", self.n
    else:
      print self.e[n-1]

parser = xml.sax.make_parser()
handler = QboxOutputHandler()
parser.setContentHandler(handler)
parser.parse(infile)
