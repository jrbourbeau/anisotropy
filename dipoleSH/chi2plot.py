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

    chi2list = ['standard', 'tibetonly', 'RI']
    lmax = 7
    chi2array = np.zeros((len(chi2list),lmax-1))
    labels = {'standard':'IceCube only', 'RI':'Combined', 'tibetonly':'Tibet Only'}

    fig = plt.figure()
    ax = plt.gca()
    redchi2dict = {}
    for chi2 in chi2list:
        redchi2list = []
        for l in range(1,lmax+1):
            # Load appropreiate SH coeff dictionary
            p = np.load('SHCoeff1e-15_1e-06/{}chi2_lmax_{}.npy'.format(chi2,l))
            p = p.item()
            redchi2list.append(p['chi2']/p['ndof'])
        plt.plot(range(1,lmax+1),redchi2list,label=labels[chi2])
        redchi2dict[chi2] = redchi2list
    print(tabulate(redchi2dict,headers="keys"))
    # print(tabulate(redchi2dict,headers="keys",tablefmt="latex"))
    plt.title(r'Reduced $\chi^2$ for 2D Fit with l_{}'.format('{max}'))
    ax.set_xlabel(r'$l_{}$'.format('{max}'))
    ax.set_ylabel(r'$\chi^2/$ndf')
    plt.legend()
    # plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    #        ncol=3, mode="expand", borderaxespad=0.)
    # plt.savefig('/home/jbourbeau/public_html/figures/2Dfit/chi2comparisons.png')
    plt.savefig('/home/jbourbeau/public_html/figures/2Dfit/chi2comparisons_optimal.png')
    ax.set_ylim(0.,15.)
    plt.savefig('/home/jbourbeau/public_html/figures/2Dfit/chi2comparisons_optimal_zoomed.png')
    # plt.savefig('/home/jbourbeau/public_html/figures/2Dfit/chi2comparisons_{}.png'.format(chi2))
