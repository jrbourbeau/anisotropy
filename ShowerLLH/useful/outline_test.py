#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('/home/jbourbeau/ShowerLLH/useful')
import geo
#data = np.array(geo.loadGeometryTxt("/net/user/santander/gamma/geometry/IT73_Outline.dat"))
#data = np.load('IT81_outline_bare.npy')
#data  = np.array(geo.loadGeometryTxt("/home/jbourbeau/ShowerLLH/useful/IC86_outline.txt"))
data = np.load('IC86_outline.npy')
labels = [format(i) for i in range(len(data))]
# print(data)
print(data.shape)
print(data.dtype)
for i in range(len(data)):
    print(format(i))
plt.scatter(data[:, 0], data[:, 1], marker = 'o', c='red', cmap = plt.get_cmap('Spectral'))

for label, x, y in zip(labels, data[:, 0], data[:, 1]):
        plt.annotate(
            label, 
            xy = (x, y), xytext = (-20, 20),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

plt.xlim(-600,600)
plt.ylim(-600,600)
plt.show()
