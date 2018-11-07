#!/usr/bin/python
# Copyright 2018 The Regents of the University of California
# This file is part of Qbox
#
# qso_plot_species.py: plot QSO potential in Gnuplot format
# using SAX incremental parsing
#
# use: qso_plot_species.py {file|URL}
import os.path
import xml.sax
import sys
import urllib2

def usage():
  print "use: ",sys.argv[0]," {file|URL}"
  sys.exit()

argc=len(sys.argv)
if ( argc != 2 ):
  usage()

input_source = sys.argv[1]

# Qbox output handler to extract and process data
class QSOSpeciesHandler(xml.sax.handler.ContentHandler):
  def __init__(self):

    self.inDescription = False
    self.inSymbol = False
    self.inAtomicNumber = False
    self.inMass = False
    self.inValenceCharge = False
    self.inMeshSpacing = False
    self.inProjector = False

    self.inNCP = False
    self.inLmax = False
    self.inLLocal = False
    self.inNquad = False
    self.inRquad = False
    self.inRadialPotential = False
    self.inRadialFunction = False

    self.inNCSLP = False
    self.inLocalPotential = False
    self.inDij = False

    self.buffer=""

  def startElement(self, name, attributes):
    if name == "symbol":
      self.symbol = ""
      self.buffer = ""
    elif name == "atomic_number":
      self.atomic_number = 0
      self.buffer = ""
    elif name == "valence_charge":
      self.valence_charge= 0
      self.buffer = ""
    elif (name == "mass"):
      self.mass = 0.0
      self.buffer = ""
    elif (name == "mesh_spacing"):
      self.mesh_spacing = 0.0
      self.buffer = ""
    elif (name == "norm_conserving_semilocal_pseudopotential"):
      self.inNCSLP = True
      print "# Norm-conserving semilocal pseudopotential"
    elif (name == "norm_conserving_pseudopotential"):
      self.inNCP = True
      print "# Norm-conserving pseudopotential"
    elif (name == "projector"):
      if self.inNCP:
        self.l = int(attributes["l"])
        self.size = int(attributes["size"])
        print "# projector l=",self.l," size=",self.size
      if self.inNCSLP:
        self.l = int(attributes["l"])
        self.i = int(attributes["i"])
        self.size = int(attributes["size"])
        print "# projector l=",self.l," i=",self.i," size=",self.size
      self.buffer = ""
    elif (name == "local_potential"):
      self.size = int(attributes["size"])
      print "# local potential, size=",self.size
      self.buffer = ""
    elif (name == "radial_potential"):
      self.buffer = ""

  def characters(self, data):
      self.buffer += data

  def endElement(self, name):
    if (name == "symbol"):
      print "# symbol:",self.buffer
      self.buffer = ""
    elif name == "atomic_number":
      print "# Z:",self.buffer
      self.buffer = ""
    elif name == "valence_charge":
      self.valence_charge=int(self.buffer)
      print "# valence charge:",self.valence_charge
      self.buffer = ""
    elif (name == "mass"):
      print "# mass:",self.buffer
      self.buffer = ""
    elif (name == "mesh_spacing"):
      self.mesh_spacing = float(self.buffer)
      print "# mesh spacing:",self.mesh_spacing
      self.buffer = ""
    elif (name == "norm_conserving_semilocal_pseudopotential"):
      self.inNCSLP = False
    elif (name == "norm_conserving_pseudopotential"):
      self.inNCP = False
    elif (name == "local_potential"):
      self.p = self.buffer.split()
      for i in range(len(self.p)):
        r = i * self.mesh_spacing
        val = float(self.p[i])
        print '%.6f'%r, '%.10e'%val
      print
      print
      self.buffer = ""
    elif (name == "projector"):
      if self.inNCSLP:
        self.p = self.buffer.split()
        for i in range(len(self.p)):
          r = i * self.mesh_spacing
          val = float(self.p[i])
          print '%.6f'%r, '%.10e'%val
        print
        print
        self.buffer = ""
    elif (name == "radial_potential"):
      print "# radial potential"
      self.p = self.buffer.split()
      for i in range(len(self.p)):
        r = i * self.mesh_spacing
        val = float(self.p[i])
        print '%.6f'%r, '%.10e'%val
      print
      print
      self.buffer = ""

parser = xml.sax.make_parser()
handler = QSOSpeciesHandler()
parser.setContentHandler(handler)
# test if input_source is a local file
# if not, process as a URL
if ( os.path.isfile(input_source) ):
  file = open(input_source)
  s = file.read(8192)
  while ( s !="" ):
    parser.feed(s)
    s = file.read(8192)
  file.close()
else:
  # attempt to open as a URL
  try:
    f = urllib2.urlopen(input_source)
    s = f.read(8192)
    while ( s !="" ):
      parser.feed(s)
      s = f.read(8192)
    f.close()
  except (ValueError,urllib2.HTTPError) as e:
    print e
    sys.exit()

parser.reset()
