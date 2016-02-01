#!/usr/bin/env python

import healpy as H
import sys
from numpy import *
from pylab import * 
from mpl_toolkits.axes_grid.parasite_axes import SubplotHost
import matplotlib as mpl

if __name__ == "__main__":

   map = H.read_map(sys.argv[1])
   frac = sum(map/max(map)) / float(len(map))

   print "Map fraction:", frac



