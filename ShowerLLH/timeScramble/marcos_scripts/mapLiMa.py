#!/usr/bin/env python

import healpy as H
import sys
from numpy import *

def LMSignificance(Non, Nback):
  alpha = 1./20.
  Noff = Nback / alpha
  sign = 1
  if Non < alpha * Noff:
     sign = -1
                  
  return sign * sqrt(2*(Non *log(((1+alpha)*Non)/(alpha * (Non + Noff))) + Noff *log(((1+alpha)*Noff)/(Non+Noff))));


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
     if map2[i] != 0 and map1[i] != 0:
        omap[i] = LMSignificance(map1[i],map2[i])

   
   H.write_map("mlima.fits", omap)


