#!/usr/bin/env python

import healpy as H
import sys
from numpy import *

if __name__ == "__main__":

   map1 = H.read_map(sys.argv[1])
   map2 = H.read_map(sys.argv[2])

   nside1 = H.get_nside(map1)
   nside2 = H.get_nside(map2)

   if nside1 != nside2: 
      print "ERROR: maps have different Nsides!"
      sys.exit()
      
   npix = H.nside2npix(nside1)
   omap = zeros(npix, dtype=double)

   for i in range(0,npix):
     omap[i] = map1[i] - map2[i]

   H.write_map("msub.fits", omap)


