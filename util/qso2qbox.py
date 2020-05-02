#!/usr/bin/python
# Convert <atomset> elements from quantum-simulation.org (QSO) format
# to Qbox input file
# use: qso2qbox.py [-last] {file|URL}
# Default: only the first <atomset> element is processed
# If using -last, only the last <atomset> element is processed

from qso import *
import os.path
import xml.sax
import sys
import urllib2
import datetime

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

# write Qbox input file
datestr=datetime.datetime.utcnow().isoformat()+'Z'
print "# converted",datestr,"from",input_source
print "set cell ",s.atoms.cell.a,s.atoms.cell.b,s.atoms.cell.c
for sp in s.atoms.species_list:
  print "species",sp.name,sp.href
for a in s.atoms.atom_list:
  print "atom",a.name,a.species,a.position[0],a.position[1],a.position[2],a.velocity[0],a.velocity[1],a.velocity[2]

