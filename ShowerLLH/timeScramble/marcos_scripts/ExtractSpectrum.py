#!/usr/bin/python

from pylab import *
import pyfits
from scipy import *
import sys

file = pyfits.open(sys.argv[1]);
cols = file[1].columns
tbdata = file[1].data
x = tbdata.field('l')
y = tbdata.field('Temperature C_l')

file = open("./spectrum.dat",'w')

for i in range(len(x)):
  #print x[i], y[i]
  line = str(x[i]) + "\t" + str(y[i]) + "\n"
  file.write(line)

file.close()
