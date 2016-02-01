#!/usr/bin/env python

import healpy as H
import sys
from numpy import *
import glob

if __name__ == "__main__":

   nside = 64
   npix = H.nside2npix(nside)
   omap = zeros(npix, dtype=double)

   for str in sys.argv[1:]:
      print "Adding file ", str
      files = glob.glob(str)
      
      for file in files:
         map = H.read_map(file)

         for i in range(0,npix):
            omap[i] += map[i]
      
   H.write_map("madd.fits", omap)


