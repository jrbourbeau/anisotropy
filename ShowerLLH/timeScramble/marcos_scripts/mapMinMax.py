#!/usr/bin/env python

import healpy as H
import sys
from numpy import *

if __name__ == "__main__":

   degree = pi / 180.

   map = H.read_map(sys.argv[1])
   nside = 64
   npix = H.nside2npix(nside)
   
   min = 1e30
   max = -1e30

   mintheta = 0
   minphi = 0
   
   maxtheta = 0
   maxphi = 0

   for i in range(0,npix):
     (theta, phi) = H.pix2ang(nside,i)
     if map[i] > max:
       max = map[i]
       maxtheta = theta
       maxphi = phi

     if map[i] < min:
       min = map[i]
       mintheta = theta
       minphi = phi


   print "Max value: ", max, " theta: ", 90 - maxtheta / degree , " phi ", maxphi / degree   
   print "Min value: ", min, " theta: ", 90 - mintheta / degree , " phi ", minphi / degree   



