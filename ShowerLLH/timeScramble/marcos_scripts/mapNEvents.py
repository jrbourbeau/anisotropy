#!/usr/bin/env python

import healpy as H
from numpy import *
from sys import *
import scipy.special as sp

if __name__ == "__main__":

   print "File: ", argv[1]
   map = H.read_map(argv[1])
   total = sum(map)
   print "Events in map: ", total
   print total/1e10, " x 10^10"



   
