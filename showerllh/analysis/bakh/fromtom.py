#!/usr/bin/env python

import re

points = {}
currentexp = ''

infile = open('spe_275_07.tab', 'r')
for line in infile:
	data = line.split()
	if len(data) == 7 and not data[5] == 'dN/dE':
		currentexp = data[5]
		points.update({currentexp:{'E':[], 'E dN/dE':[], 'E2.75 dN/dE':[], 'error':[]}})
		data = data[:4]
	if len(data) == 4:
		points[currentexp]['E'].append(float(data[0]))
		points[currentexp]['E dN/dE'].append(float(data[1]))
		points[currentexp]['E2.75 dN/dE'].append(float(data[2]))
		points[currentexp]['error'].append(float(data[3]))


