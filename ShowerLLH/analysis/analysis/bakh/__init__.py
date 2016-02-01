#!/usr/bin/env python

import re
import os
import numpy as np

points = {}
def dloge2de(e,dloge):
	return dloge*np.power(e,-1)/np.log(10.0)
def weight(e,dnde):
	return np.power(e,2.7)*dnde
def logweight(e,dnde):
	return np.log10(weight(e,dnde))
def logweightederror(e,dnde,low,high):
	low = logweight(e,dnde-low)
	high = logweight(e,dnde+high)
	return logweight(e,dnde)-low, high-logweight(e,dnde)
def geterrorbars(exp):
	data = points[exp]
	return logweightederror(np.array(data['E']),np.array(data['dN/dE']),np.array(data['error_low']),np.array(data['error_high']))

datadir = '/home/jeisch/compfit/others/spectrum'
for file in os.listdir(datadir):
	infile = open(os.path.join(datadir,file), 'r')
	name = file.split('.')[0]
	points.update({name:{'E':[], 'dN/dE':[], 'E2.7 dN/dE':[], 'error_high':[], 'error_low':[]}})
	for line in infile:
		if line[0] == '#':
			continue
		data = line.split()
		if not len(data) in (2,4):
			continue
		data = [float(f) for f in data]
		if name[:4] == "IT73":
			E = np.power(10.0,data[0])
			points[name]['E'].append(E)
			dnde = dloge2de(E,data[1])
			points[name]['dN/dE'].append(dnde)
			points[name]['E2.7 dN/dE'].append(weight(E,dnde))
		else:
			points[name]['E'].append(data[0])
			points[name]['dN/dE'].append(data[1])
			points[name]['E2.7 dN/dE'].append(weight(data[0],data[1]))
		if len(data) < 4:
			continue
		points[name]['error_high'].append(data[2])
		points[name]['error_low'].append(data[3])



lnadir = '/home/jeisch/compfit/others/lna/'
lna = {}
for file in os.listdir(lnadir):
	infile = open(os.path.join(lnadir,file), 'r')

	name = file.split('.')[0]
	lna.update({name:{'E':[], 'lnA':[], 'error_high':[], 'error_low':[]}})
	for line in infile:
		if line[0] == '#':
			continue
		data = line.split()
		if not len(data) in (2,4):
			continue
		data = [float(f) for f in data]
		lna[name]['E'].append(data[0])
		lna[name]['lnA'].append(data[1])
		if len(data) < 4:
			continue
		lna[name]['error_high'].append(data[2])
		lna[name]['error_low'].append(data[3])


