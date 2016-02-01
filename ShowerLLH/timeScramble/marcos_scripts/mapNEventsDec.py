#!/usr/bin/env python

import healpy as H
from numpy import *
from sys import *
import scipy.special as sp

if __name__ == "__main__":

   if len(argv) < 4: 
      print "usage: mapNEventsDec.py [mindec] [maxdec] [FITS file]"
      sys.exit()

   
   degree = pi / 180.
   thetai = (90 - float(argv[2])) * degree
   thetaf = (90 - float(argv[1])) * degree

   print "low edge: ", argv[1], " hi edge: ", argv[2]
   
   print "File: ", argv[3]
   map = H.read_map(argv[3])

   nside = H.get_nside(map)
   npix = H.nside2npix(nside)
   sum = 0

   for i in range(0,npix):
     (theta, phi) = H.pix2ang(nside,i)
     if theta >= thetai and theta <= thetaf:
        sum += map[i]
   
   print "Events in map: ", sum
   print sum/1e10, " x 10^10"



   
