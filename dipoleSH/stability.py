#!/usr/bin/env python

#=========================================================================
# File Name     : chi2plot.py
# Description   : Plot chi2
# Creation Date : 03-18-2016
# Last Modified : Fri 18 Mar 2016 12:24:44 PM CDT
# Created By    : James Bourbeau
#=========================================================================

import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate

if __name__ == "__main__":

    # chi2list = ['standard','RI']
    chi2list = ['standard', 'tibetonly', 'RI']
    lmax = 4
    labels = {'standard':'IceCube only', 'RI':'Combined', 'tibetonly':'Tibet Only'}
    fig = plt.figure()
    ax = plt.gca()
    datadict = {}
    init = 1e-15
    step = 1e-06
    for chi2 in chi2list:
        almdict = {}
        for l in range(1,lmax+1):
            # Load appropreiate SH coeff dictionary
            # p = np.load('SHcoeff0/{}chi2_lmax_{}.npy'.format(chi2,l))
            p = np.load('SHCoeff/SHCoeff{}_{}/{}chi2_lmax_{}.npy'.format(init,step,chi2,l))
            p = p.item()
            almdict[l] = [p['Y(1,-1)'],p['dY(1,-1)'],p['Y(1,1)'],p['dY(1,1)']]
        # plt.plot(range(1,lmax+1),redchi2list,label=labels[chi2])
        datadict[chi2] = almdict
    for chi2 in chi2list:
        x = np.array([datadict[chi2][i][0] for i in range(1,lmax+1)])
        dx = np.array([datadict[chi2][i][1] for i in range(1,lmax+1)])
        y = np.array([datadict[chi2][i][2] for i in range(1,lmax+1)])
        dy = np.array([datadict[chi2][i][3] for i in range(1,lmax+1)])
        plt.errorbar(10**4*x,10**4*y,xerr=10**4*dx,yerr=10**4*dy,label=labels[chi2],linestyle='-.',marker='.')
    # print(tabulate(redchi2dict,headers="keys"))
    plt.plot(0,0,marker='.',color='black')
    # plt.plot(0,0,marker='*',color='black')
    ax.set_xlabel(r'$a_{}$'.format('{1,-1}'))
    ax.set_ylabel(r'$a_{}$'.format('{1,1}'))
    scale = 40.
    ax.set_ylim(-scale,scale)
    ax.set_xlim(-scale,scale)
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=3, mode="expand", borderaxespad=0.)
    plt.savefig('/home/jbourbeau/public_html/figures/2Dfit/almmap.png')
