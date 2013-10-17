#!/usr/bin/python
# qbox_dos.py: extract dos from Qbox output
# generate dos plot in gnuplot format
# use: qbox_dos.py emin emax width file.r

import xml.sax
import sys
import math

if len(sys.argv) != 5:
  print "use: ",sys.argv[0]," emin emax width file.r"
  sys.exit()

emin = float(sys.argv[1])
emax = float(sys.argv[2])
width = float(sys.argv[3])
infile = sys.argv[4]

ndos = 501
de = (emax - emin)/(ndos-1)

# normalized gaussian distribution in one dimension
# f(x) = 1/(sqrt(pi)*width) * exp(-(x/width)^2 )
def gauss(x, width):
  return (1.0/(math.sqrt(math.pi)*width)) * math.exp(-(x/width)**2)
  
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
      self.compute_dos()

  def compute_dos(self):
    dos = [0] * ndos
    self.e = self.buffer.split()
    for i in range(len(self.e)):
      for j in range(ndos):
        ej = emin + j * de
        dos[j] += gauss(float(self.e[i])-ej, width )
    if self.iter > 1:
      print
      print
    print "# ",infile, " iter=", self.iter," spin=",self.spin, \
          " kpoint=", self.kpoint," width=",width
    for j in range(ndos):
      ej = emin + j * de
      print ej, dos[j]

parser = xml.sax.make_parser()
handler = QboxOutputHandler()
parser.setContentHandler(handler)
parser.parse(infile)
