#!/usr/bin/env python
from __future__ import division
from matplotlib.patches import Patch
from pylab import *
import numpy as np
import healpy as H

def GetIntValue(map,theta,phi):
    return H.get_interp_val(map,theta,phi)

degree = pi / 180.

N = 300
phi = np.linspace(0, 360 * degree, 2*N)
theta = np.linspace(0, 180 * degree, N)

Phi, Theta = meshgrid(phi, theta)

map = H.read_map(sys.argv[1])
nside1 = H.get_nside(map)
npix = H.nside2npix(nside1)
map = map 

palette = cm.jet
#mask = palette.set_under('g', 1.0)
mapmask = np.ma.masked_where(map == 0, map)

fig = figure(figsize=(12,6))
value = GetIntValue(mapmask,Theta, 360 * degree - Phi)
ax = subplot(111)
im = imshow(value, cmap=cm.jet,extent=[360,0,-90,90])
#im.set_interpolation('bilinear')
colorbar(orientation="horizontal",shrink=0.5,fraction=0.05,pad=0.1)
grid()
clim(-1e-3,2e-3)

show()
