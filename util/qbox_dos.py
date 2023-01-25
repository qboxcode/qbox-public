#!/usr/bin/env python3
# qbox_dos.py: extract electronic DOS from Qbox output
# generate DOS plot in gnuplot format
# use: qbox_dos.py [-last] [mu] emin emax [mu] width file.r
# mu: (optional) chemical potential [eV]
#     If omitted, an attempt is made to read the element <mu> from file.r
#     The value is 0.0 if <mu> is not found in file.r
# emin, emax: bounds of plot in [eV] relative to chemical potential mu
# width: gaussian broadening width [eV]
# the DOS is accumulated separately for each spin
# With the -last option, only the last <eigenset> is used to compute the DOS

import xml.sax
import sys
import math

argc=len(sys.argv)
if (argc < 5) or (argc > 7) :
  print("use: ",sys.argv[0]," [-last] [mu] emin emax width file.r")
  sys.exit()

# default chemical potential mu=0.0
mu=0.0

iarg = 1
lastonly = False
# check for -last option
if (sys.argv[iarg] == "-last") :
  lastonly = True
  iarg += 1
  if (argc < 6):
    print("use: ",sys.argv[0]," [-last] [mu] emin emax width file.r")
    sys.exit()

# check for mu argument
if ((lastonly and (argc == 7)) or ((not lastonly) and (argc == 6))):
  mu = float(sys.argv[iarg])
  iarg += 1
  if (argc < 6):
    print("use: ",sys.argv[0]," [-last] [mu] emin emax width file.r")
    sys.exit()

emin = float(sys.argv[iarg])
iarg += 1
emax = float(sys.argv[iarg])
iarg += 1
width = float(sys.argv[iarg])
iarg += 1
infile = sys.argv[iarg]

ndos = 501
de = (emax - emin)/(ndos-1)

# normalized gaussian distribution in one dimension
# f(x) = 1/(sqrt(pi)*width) * exp(-(x/width)^2 )
def gauss(x, width):
  return (1.0/(math.sqrt(math.pi)*width)) * math.exp(-(x/width)**2)

# Qbox output handler to extract and process data
class QboxOutputHandler(xml.sax.handler.ContentHandler):
  def __init__(self):
    self.nspin = 1
    self.readData = 0
    self.dos_up = [0] * ndos
    self.dos_dn = [0] * ndos
    self.mu = mu

  def startElement(self, name, attributes):
    if (name == "eigenset") and (lastonly):
      self.dos_up = [0] * ndos
      self.dos_dn = [0] * ndos
    if name == "eigenvalues":
      self.n = attributes["n"]
      self.spin = int(attributes["spin"])
      self.kpoint = attributes["kpoint"]
      self.weight = float(attributes["weight"])
      self.readData = 1
      self.buffer = ""
      if self.spin == 1:
        self.nspin = 2
    if name == "mu":
      self.readData = 1
      self.buffer=""

  def characters(self, data):
    if self.readData:
      self.buffer += data

  def endElement(self, name):
    if name == "eigenvalues":
      self.readData = 0
      self.accumulate_dos()
    if name == "mu":
      self.readData = 0
      self.mu = float(self.buffer)

  def accumulate_dos(self):
    self.e = self.buffer.split()
    if self.spin == 0:
      for i in range(len(self.e)):
        for j in range(ndos):
          ej = emin + j * de
          self.dos_up[j] += gauss(float(self.e[i])-self.mu-ej, width ) * self.weight
    if self.spin == 1:
      for i in range(len(self.e)):
        for j in range(ndos):
          ej = emin + j * de
          self.dos_dn[j] += gauss(float(self.e[i])-self.mu-ej, width ) * self.weight

  def print_dos(self):
    print("# ",infile," mu=",self.mu," spin=0 width=",width)
    for j in range(ndos):
      ej = emin + j * de
      print(ej, self.dos_up[j])
    if self.nspin == 2:
      print()
      print()
      print("# ",infile," mu=",mu," spin=1 width=",width)
      for j in range(ndos):
        ej = emin + j * de
        print(ej, self.dos_dn[j])

parser = xml.sax.make_parser()
handler = QboxOutputHandler()
parser.setContentHandler(handler)
parser.parse(infile)
handler.print_dos()
