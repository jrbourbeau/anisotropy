#!/usr/bin/env python

import healpy as H
import sys
from numpy import *
from pylab import * 
from mpl_toolkits.axes_grid.parasite_axes import SubplotHost
import matplotlib as mpl
import matplotlib.mlab as mlab

if __name__ == "__main__":

   fig = figure()
   ax = fig.add_subplot(111)

   imap = H.read_map(sys.argv[1])
   map = [x for x in imap if x != 0]
   n, bins, patches = ax.hist(map, 50, normed=1, facecolor='green', alpha=0.75)

   show() 


