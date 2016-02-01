#!/usr/bin/env python

import healpy as H
import sys
from numpy import *

if __name__ == "__main__":

   print "./cmd [rima] [datamap]"

   rimap = H.read_map(sys.argv[1])
   datamap = H.read_map(sys.argv[2])

   nside1 = H.get_nside(rimap)

   npix = H.nside2npix(nside1)
   omap = zeros(npix, dtype=double)
   
   for i in range(0,npix):
       omap[i] = datamap[i] / (rimap[i] + 1)

   H.write_map("mbkgri.fits", omap)
   

