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

    initlist = [0.,1e-15,1e-12,1e-9,1e-6]
    steplist = [1e-6,1e-9,1e-12]

    minchi2dict = {}
    for chi2 in chi2list:
        fig = plt.figure()
        ax = plt.gca()
        # redchi2dict = {}
        if chi2!='tibetonly':
            minchi2dict[chi2] = [1000,1000,1000]
        for init in initlist:
            for step in steplist:
            # for chi2 in chi2list:
                redchi2list = []
                for l in range(1,lmax+1):
                    # Load appropreiate SH coeff dictionary
                    p = np.load('SHCoeff/SHCoeff{}_{}/{}chi2_lmax_{}.npy'.format(init,step,chi2,l))
                    p = p.item()
                    redchi2list.append(p['chi2']/p['ndof'])
                # plt.plot(range(1,lmax+1),redchi2list,label=labels[chi2])
                if max(redchi2list) >= 50 and chi2!='tibetonly':
                    print('\nBad params:')
                    print('\tchi2 = {}'.format(chi2))
                    print('\tinit = {}'.format(init))
                    print('\tstep = {}\n'.format(step))
                if chi2!='tibetonly' and min(redchi2list)<=minchi2dict[chi2][0]:
                    minchi2dict[chi2][0] = min(redchi2list)
                    minchi2dict[chi2][1] = init
                    minchi2dict[chi2][2] = step
                plt.plot(range(1,lmax+1),redchi2list,label='{}/{}'.format(init,step))
                # plt.plot(range(1,lmax+1),redchi2list,label='{}/{}/{}'.format(labels[chi2],init,step))
                # redchi2dict[chi2] = redchi2list
        # print(tabulate(redchi2dict,headers="keys",tablefmt="latex"))
        # plt.title(r'Reduced $\chi^2$ for 2D Fit with l_{}'.format('{max}'))
        ax.set_xlabel(r'$l_{}$'.format('{max}'))
        ax.set_ylabel(r'$\chi^2/$ndf')
        # plt.legend()
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=3, mode="expand", borderaxespad=0.)
        # plt.savefig('/home/jbourbeau/public_html/figures/2Dfit/chi2comparisons.png')
        # ax.set_ylim(0.,31.)
        # plt.savefig('/home/jbourbeau/public_html/figures/2Dfit/chi2comparisons_zoomed.png')
        # plt.savefig('/home/jbourbeau/public_html/figures/2Dfit/chi2comparisons_{}.png'.format(chi2))
    print('minchi2dict = {}'.format(minchi2dict))
    # plt.show()
