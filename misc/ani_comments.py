#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import argparse
import os, re, sys
import colormaps as cmaps
from useful import getMids

if __name__ == "__main__":
    
    public_path = '/home/jbourbeau/public_html/figures/'
    parser = argparse.ArgumentParser(description='Calculate anisotropy')
    parser.add_argument('--nbins', dest='num_bins', type=int, default=6, help='Number of points to be used')
    parser.add_argument('--lowerb', dest='lower_bound', default=False, action='store_true', help='Use extended lower bound')
    parser.add_argument('--out', dest='out', default=False, action='store_true', help='Save plots')
    parser.add_argument('--weighted', dest='weighted', default=False, action='store_true', help='Use weighted histograms over all IC configs')
    parser.add_argument('--config', dest='config', default=None, help='Detector configuration')

    args = parser.parse_args()
    
    # Basic setup
    xmin = 2.75
    xmax = 9.01
    step_size = 0.05
    number_steps = (xmax-xmin)/step_size
    
    energy_sim_bins = np.arange(xmin, xmax, step_size)
    x = getMids(energy_sim_bins)
    rel_freq = [0.0]*9
    
    config_events = {'IC59':3.579e+10, 'IC79':4.131e+10, 'IC86':5.906e+10, 'IC86-II':5.630e+10, \
            'IC86-III':6.214e+10, 'IC86-IV':6.327e+10}
    total = 2.5463e+11
    
    # Read in relative frequency vs. energy histograms
    fig, ax = plt.subplots()
    if args.weighted:
        args.config = 'weighted'
        configuration = ['IC59','IC79','IC86','IC86-II','IC86-III']
        for Config in configuration:
            f = np.load('/home/fmcnally/anisotropy/icesim/{}_hists.npy'.format(Config))
            for i, y in enumerate(f):
                #print('i = {}'.format(i))
                total_num_events = float(y.sum())
                rel_freq[i] += np.array((y/total_num_events)*(config_events[Config]/total))
        for j in range(len(rel_freq)):
            #print('weighted_{} = {}'.format(j+1,rel_freq[j]))
            ax.step(x, rel_freq[j], label='{}'.format(i+1))
    
    else:
        if args.config == None:
            raise SystemExit('\nError: Either a detector configuration or the wegithed option must be specified\n')
        
        f = np.load('/home/fmcnally/anisotropy/icesim/{}_hists.npy'.format(args.config))
        for i, y in enumerate(f):
            total_num_events = float(y.sum())
            ax.step(x, y/total_num_events, label=i+1)
            rel_freq[i] = y/total_num_events
            #print('{}_{} = {}'.format(args.config,i+1,y/total_num_events))
    
    ax.set_xlim(2.75, 8.5)
    tPars = {'fontsize':16}
    ax.set_xlabel(r'True Energy ($\log_{10}(E/\mathrm{GeV})$)', **tPars)
    ax.set_ylabel('Fraction of Events', **tPars)
    plt.legend(loc='upper right')
    plt.savefig(public_path+'rel_freq_{}'.format(args.config), dpi=300, bbox_inches='tight')
    
    ## Quickly calculate the median
    #values = []
    #for i in range(len(x)):
    #    print('[x[i]]*rel_freq[0][i] = {}'.format([x[i]]*rel_freq[0][i]))
    #    values.append([x[i]]*rel_freq[0][i])
    #print('median = {}'.format(np.median([item for sublist in values for item in sublist])))
    
    # Calculate desired anisotropy function bin widths 
    # Derived from the median energy points in Figure 8 		
    energy_bin_median = [4.12, 4.38, 4.58, 4.85, 5.12, 5.38, 5.77, 6.13, 6.73]
    ani_bins = [0.0]*len(energy_bin_median)
    for i in range(len(energy_bin_median)):
        if i == 0:
            bin_upper = (energy_bin_median[i+1]-energy_bin_median[i])/2.0
            bin_lower = bin_upper
            if not args.lower_bound:
                ani_bins[i] = ([energy_bin_median[i]-bin_lower,energy_bin_median[i]+bin_upper])
            else:
                ani_bins[i] = ([3.0,energy_bin_median[i]+bin_upper])
        if i == len(energy_bin_median)-1:
            bin_lower = (energy_bin_median[i]-energy_bin_median[i-1])/2.0
            bin_upper = bin_lower
            ani_bins[i] = ([energy_bin_median[i]-bin_lower,energy_bin_median[i]+bin_upper])
        if i != 0 and i != len(energy_bin_median)-1:
            bin_lower = (energy_bin_median[i]-energy_bin_median[i-1])/2.0
            bin_upper = (energy_bin_median[i+1]-energy_bin_median[i])/2.0
            ani_bins[i] = ([energy_bin_median[i]-bin_lower,energy_bin_median[i]+bin_upper])
    
    # Calculate the N matrix elements
    def get_matrix_element(i,j):
        bin_min = ani_bins[j][0]
        bin_max = ani_bins[j][1]
        element_ij = 0.0
        for x in range(len(energy_sim_bins)):
            if energy_sim_bins[x] >= bin_min and energy_sim_bins[x] < bin_max:
                element_ij += rel_freq[i][x]
        return element_ij
    
    N = []
    for i in range(args.num_bins):
        row = []
        for j in range(args.num_bins):
            row.append(get_matrix_element(i,j))
        N.append(row)
    
    '''plt.figure(2)
    mat = plt.matshow(np.matrix(N).getI(),fignum=2, cmap=cmaps.viridis)
    #mat = plt.matshow(np.matrix(N).getI(),fignum=2, cmap=plt.get_cmap('coolwarm'))
    plt.colorbar(mat, orientation='vertical')
    plt.figure(4)
    im = plt.matshow(np.matrix(N), fignum=4, cmap=cmaps.viridis)
    #im = plt.matshow(np.matrix(N), fignum=4, cmap=plt.get_cmap('coolwarm'))
    plt.colorbar(im, orientation='vertical')'''

    anisotropy_data = [7.32170723e-04, 6.96619874e-04, 4.82783053e-04, 2.85341450e-04, 1.60042336e-04, 7.31462195e-05, 2.00920356e-04, 5.00927901e-04, 6.52328253e-04]
    data_used = []
    for i in range(args.num_bins): 
        data_used.append(anisotropy_data[i])
    sol = np.linalg.solve(np.array(N),np.array(data_used))
    #print(sol)

    def d_I(E):
        for i in range(args.num_bins):
            lower = ani_bins[i][0]
            upper = ani_bins[i][1]
            #if E < ani_bins[0][0] or E > ani_bins[len(ani_bins)-1][1]:
            #	return 'Not in range.'
            if E >= lower and E < upper:
                return sol[i]
    
    plt.figure(3)
    for i in range(args.num_bins):
        lower = ani_bins[i][0]
        upper = ani_bins[i][1]
        xnew = np.arange(lower, upper, 0.01);
        plt.plot(xnew, [d_I(x) for x in xnew], marker='None', linestyle='-', linewidth=4.0, color='r') 
    plt.xlabel(r'True Energy [$\log_{10}$(E/GeV)]', **tPars)
    plt.ylabel(r'$\delta I(E)$', **tPars)
    title = public_path+'function_{}_nbins{}'.format(args.config,args.num_bins)
    if args.lower_bound:
        title += '_lowerb'
    plt.savefig(title, dpi=300, bbox_inches='tight')
    plt.show()
