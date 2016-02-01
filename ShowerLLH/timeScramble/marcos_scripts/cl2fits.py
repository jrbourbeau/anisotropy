#!/usr/bin/env python

from pylab import *
import numpy as n
import pyfits
import sys

clfile = n.genfromtxt(sys.argv[1])
l = clfile[:,0]
Cl = clfile[:,1]

newl = pyfits.Column(name='l',format='J', array=l)
newCl = pyfits.Column(name='Temperature C_l', format='D',array=Cl)

#ascii table
tbhdu=pyfits.new_table([newCl],tbtype='TableHDU')
tbhdu.writeto('spectrum.fits')

