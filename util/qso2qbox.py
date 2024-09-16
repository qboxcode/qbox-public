#!/usr/bin/env python3
# Convert <atomset> elements from quantum-simulation.org (QSO) format
# to Qbox input file
# use: qso2qbox.py [-last] {file|URL}
# Default: only the first <atomset> element is processed
# If using -last, only the last <atomset> element is processed

from qso import *
import os.path
import xml.sax
import sys
import urllib
import datetime

def usage():
  print ("use: ",sys.argv[0]," [-last] {file|URL}")
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
    f = urllib.request.urlopen(input_source)
    str = f.read(8192)
    while ( str !="" and not (first_only and handler.done_first) ):
      parser.feed(str)
      str = f.read(8192)
    f.close()
  except (ValueError,urllib.error.HTTPError) as e:
    print (e)
    sys.exit()

parser.reset()

# write Qbox input file
datestr=datetime.datetime.utcnow().isoformat()+'Z'
print ("# converted",datestr,"from",input_source)
av = s.atoms.cell.a
bv = s.atoms.cell.b
cv = s.atoms.cell.c
print ("set cell ",
 '%.8f'%(av[0]), '%.8f'%(av[1]), '%.8f'%(av[2]),
 '%.8f'%(bv[0]), '%.8f'%(bv[1]), '%.8f'%(bv[2]),
 '%.8f'%(cv[0]), '%.8f'%(cv[1]), '%.8f'%(cv[2]))
for sp in s.atoms.species_list:
  print ("species",sp.name,sp.href)
for a in s.atoms.atom_list:
  print ("atom",a.name,a.species,'%14.8f'%a.position[0],'%14.8f'%a.position[1],'%14.8f'%a.position[2],'%14.8f'%a.velocity[0],'%14.8f'%a.velocity[1],'%14.8f'%a.velocity[2])
