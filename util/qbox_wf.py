#!/usr/bin/env python3
# extract a Qbox wave function from a base64 encoded restart file
# save wf in a numpy 3d array
# save wf in an npy file
import sys
import numpy as np
from qso import getwf

if len(sys.argv) != 3:
  print ("use: ",sys.argv[0]," file.xml n")
  sys.exit()

filename=sys.argv[1]
n=int(sys.argv[2])
w=getwf(filename,n)
np.save("wf.npy",w)
