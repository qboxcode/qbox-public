#!/usr/bin/python
# Convert from quantum-simulation.org format to QuantumEspresso input format
# use: env ECUT=ecut qso2qe.py [-last] {file|URL}

from qso import *
import sys
import os.path
import urllib2

# Get Ecut from environment variable, default=25
ecut = os.getenv("ECUT",25.0)

def usage():
  print "use: ",sys.argv[0]," [-last] {file|URL}"
  sys.exit()

argc=len(sys.argv)
if ( argc < 2 or argc > 3 ):
  usage()

# check if option "-last" is used
# "-last" option: process all <atomset> and return only the last
# default: extract first atomset only
first_only = True
input_source = sys.argv[1]
if ( sys.argv[1] == "-last" ):
  if ( argc != 3 ):
    usage()
  first_only = False
  input_source = sys.argv[2]

s = Sample()
parser = xml.sax.make_parser()
handler = QSOAtomSetHandler(s)
parser.setContentHandler(handler)
# test if input_source is a local file
# if not, process as a URL
if ( os.path.isfile(input_source) ):
  file = open(input_source)
  str = file.read(8192)
  while ( str !="" and not (first_only and handler.done_first) ):
    parser.feed(str)
    str = file.read(8192)
  file.close()
else:
  # attempt to open as a URL
  try:
    f = urllib2.urlopen(input_source)
    str = f.read(8192)
    while ( str !="" and not (first_only and handler.done_first) ):
      parser.feed(str)
      str = f.read(8192)
    f.close()
  except (ValueError,urllib2.HTTPError) as e:
    print e
    sys.exit()

parser.reset()

# write QE input file

print "&control"
print "  calculation = 'scf'"
print "  pseudo_dir = './'"
print "/"
print "&system"
print "  ibrav=0"
print "  nat=",len(s.atoms.atom_list),", ntyp=",len(s.atoms.species_list),","
print "  ecutwfc=",ecut
print "/"
print "&electrons"
print "/"
print "ATOMIC_SPECIES"
for sp in s.atoms.species_list:
  print sp.name,sp.mass,sp.href
print "ATOMIC_POSITIONS {bohr}"
for a in s.atoms.atom_list:
  print a.species,a.position[0],a.position[1],a.position[2]
print "CELL_PARAMETERS bohr"
print s.atoms.cell.a
print s.atoms.cell.b
print s.atoms.cell.c
print "K_POINTS {gamma}"
