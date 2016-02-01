#!/usr/bin/env python

import numpy as np

e = np.power(10, 5.5+np.arange(13)/5.0)

h = np.array([2.60e6, 4.32e6, 2.86e6, 1.14e6, 503145, 177977, 118716, 40223, 9982, 2545, 2166, 65, 218])
he = np.array([2.93e6, 1.82e6, 1.49e6, 0.82e6, 557799, 250468, 97561, 28340, 14347, 4908, 2342, 262, 208])
o = np.array([2.58e6, 0.85e6, 0.81e6, 0.51e6, 280927, 96900, 25476, 15600, 7109, 4309, 1200, 649, 69])
si = np.array([1.66e6, 1.36e6, 0.87e6, 0.59e6, 183970, 85332, 36698, 18283, 4821, 3897, 705, 679, 879])
fe = np.array([2.33e6, 1.59e6, 1.31e6, 0.25e6, 63928.4, 28479, 6179, 5000, 506, 943, 828, 777, 254])
total = np.array([12.1e6, 9.94e6, 7.34e6, 3.32e6, 1.59e6, 639156, 284630, 107445, 36765, 16653, 6496, 2430, 1629])

lna = h/total*np.log(1) + he/total*np.log(4)+o/total*np.log(16)+si/total*np.log(28)+fe/total*np.log(56)

f = open('lna/chen.dat', 'w')

for i in range(len(e)):
	f.write('{0}\t{1}\n'.format(e[i], lna[i]))

f.close()


