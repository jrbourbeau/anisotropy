#!/usr/bin/env python

import healpy as H
import sys
from numpy import *

if __name__ == "__main__":

   degree = pi / 180.

   ra = [ 122.4, 263.0, 201.6, 332.4, 217.7, 77.6, 308.2, 166.5 ]
   dec = [ -47.4, -44.1, -37.0, -70.0, -70.0, -31.9, -34.5, -37.2 ]
   num = [ '1', '2', '3', '4', '5', '6', '7', '8']
   scale = [22, 13, 11, 12, 12, 13, 20, 12]
   sigma = [7.0, 6.7, 6.3, 6.2, -6.4, -6.1, -6.1, -6.0]

   for i in range(0,len(ra)):
     filename = "ic79_dec25_signal_%ddeg.fits" % scale[i]

     theta = (90 - dec[i]) * degree
     phi = ra[i] * degree

     map1 = H.read_map(filename)
     nside = H.get_nside(map1)

     pix = H.ang2pix(nside, theta, phi)
     sig = '%.1f' % map1[pix]

     print "#", i+1, " RA: ", ra[i], " Dec: ", dec[i], " Scale: ", scale[i], "Significance: ", sig, " was: ", sigma[i]


