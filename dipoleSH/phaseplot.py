#!/usr/bin/env python

#=========================================================================
# File Name     : phaseplot.py
# Description   : Plot phase
# Creation Date : 03-18-2016
# Last Modified : Fri 18 Mar 2016 01:48:49 PM CDT
# Created By    : James Bourbeau
#=========================================================================

import numpy as np
import matplotlib.pyplot as plt
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Calculate phase from 2D fit map')
    parser.add_argument('--init', dest='init', type=float, \
        help='Number of RA bins to be used in projection')
    parser.add_argument('--step', dest='step', type=float, \
        help='Number of RA bins to be used in projection')

    args = parser.parse_args()
    opts = vars(args).copy()

    init = opts['init']
    step = opts['step']

    chi2list = ['standard', 'tibetonly', 'RI']
    labels = {'standard':'IceCube only', 'RI':'Combined', 'tibetonly':'Tibet only'}
    fig = plt.figure()
    ax = plt.gca()
    for chi2 in chi2list:
        phaselist = []
        errlist = []
        for l in range(1,8):
            # Load appropreiate SH coeff dictionary
            p = np.load('SHCoeff{}_{}/{}chi2_lmax_{}.npy'.format(init,step,chi2,l))
            p = p.item()
            phaselist.append(np.arctan(p['Y(1,1)']/p['Y(1,-1)'])*180./np.pi)

            prefactor = 1./(1+(p['Y(1,1)']/p['Y(1,-1)'])**2)
            first = (p['dY(1,1)']/p['Y(1,-1)'])**2
            second = ((p['dY(1,-1)']*p['Y(1,1)'])/(p['Y(1,-1)']**2))**2
            error = prefactor*np.sqrt(first+second)*180./np.pi
            errlist.append(error)
            # print('\nl/chi2 = {}/{}'.format(l,chi2))
            # print('ierr = {}\n'.format(p['ierr']))
            # if abs(error) >= 50:
            #     print('\nl/chi2 = {}/{}'.format(l,chi2))
            #     print('prefactor = {}'.format(prefactor))
            #     print('first = {}'.format(first))
            #     print('second = {}'.format(second))
            #     print('error = {}'.format(error))
            #     print('Y(1,1) = {}'.format(p['Y(1,1)']))
            #     print('Y(1,-1) = {}'.format(p['Y(1,-1)']))
            #     print('dY(1,1) = {}'.format(p['dY(1,1)']))
            #     print('dY(1,-1) = {}'.format(p['dY(1,-1)']))
            #     print('ierr = {}\n'.format(p['ierr']))
            x = np.arange(1,8)
            if chi2=='tibetonly':
                x = x-0.05
                x = x[:3]
                phaselist = phaselist[:3]
                errlist = errlist[:3]
            if chi2=='RI':
                x = x+0.05
        plt.errorbar(x,phaselist,yerr=errlist,marker='.',linestyle='-.',label=labels[chi2]+r' $\chi^2$',alpha=0.7)
        # plt.plot(range(1,8),phaselist,label=labels[chi2]+r' $\chi^2$')
    # plt.title(r'Phase from 2D fit maps')
    ax.set_xlabel(r'$l_{}$'.format('{max}'))
    ax.set_ylabel(r'Dipole phase $[^{}]$'.format('{\circ}'))
    # plt.legend()
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=3, mode="expand", borderaxespad=0.)
    ax.set_ylim(-100,100)
    ax.set_xlim(0.75,7.25)
    plt.savefig('/home/jbourbeau/public_html/figures/2Dfit/dipole-phase_SHCoeff{}_{}.png'.format(init,step))
    # plt.show()
