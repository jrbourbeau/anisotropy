#!/usr/bin/env python

import healpy as H
import sys
from numpy import *
from math import *

if __name__ == "__main__":

   map1 = H.read_map(sys.argv[1])

   nside1 = H.get_nside(map1)

   npix = H.nside2npix(nside1)
   omap = zeros(npix, dtype=double)

   for i in range(0,npix):
     sign = map1[i] / fabs(map1[i])
     omap[i] = pow(fabs(map1[i]), float(sys.argv[2])) * sign


     
     #,  * map1[i] * sign

   ofile = "msqr.fits"
   H.write_map(ofile, omap)


