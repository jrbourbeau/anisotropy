#!/usr/bins/env python

import numpy as np

pf = '/data/user/zgriffith/ShowerLLH/resources/'
lt = np.load(pf+'LLHTables.npy')
b  = np.load(pf+'ShowerLLH_bins.npy')

llhtables   = lt.item()
bins        = b.item()
new_bins    = {}
new_bins[0] = ['E',bins['Ebins']]
new_bins[1] = ['Z',bins['Zbins']]
new_bins[2] = ['S',bins['Sbins']]
new_bins[3] = ['D',bins['Dbins']]
new_bins[4] = ['C',bins['Cbins']]

new_tables           = {}
new_tables['proton'] = llhtables['P']
new_tables['iron' ]  = llhtables['Fe']

new_file              = {}
new_file['bins']      = new_bins
new_file['llhtables'] = new_tables

np.save(pf+'new_LLHTables.npy', new_file)

