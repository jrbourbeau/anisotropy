#!/usr/bin/env python

#==============================================================================
# File Name     : SHcomparison.py
# Description   : Compare the old (hard-coded) and new
#                 (scipy) versions of generating SH
# Creation Date : 03-08-2016
# Last Modified : Tue 08 Mar 2016 01:19:22 PM CST
# Created By    : James Bourbeau
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib import rc
import os, glob

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

if __name__ == "__main__":
    oldparamfiles = glob.glob('/home/jbourbeau/anisotropy/dev/SHcomparison/old*.npy')
    oldparamfiles.sort()
    newparamfiles = glob.glob('/home/jbourbeau/anisotropy/dev/SHcomparison/new*.npy')
    newparamfiles.sort()

    oldparamdict = {}
    for i in range(len(oldparamfiles)):
        l = oldparamfiles[i][-5:-4]
        oldparamdict['{}'.format(l)] = np.load(oldparamfiles[i]).item()
    newparamdict = {}
    for i in range(len(newparamfiles)):
        l = newparamfiles[i][-5:-4]
        newparamdict['{}'.format(l)] = np.load(newparamfiles[i]).item()

    # Useful breakdown of l, m values to be used
    l=4
    nsph = sum([2*l_i+1 for l_i in range(l+1)])
    lvals = [[l_i]*(2*l_i+1) for l_i in range(l+1)]
    mvals = [[m for m in range(-l_i, l_i+1)] for l_i in range(l+1)]
    lvals = [item for sublist in lvals for item in sublist]
    mvals = [item for sublist in mvals for item in sublist]
    fitparams = ['Y(%i,%i)' % (lvals[i], mvals[i]) for i in range(nsph)]
    fiterrparams = ['dY(%i,%i)' % (lvals[i], mvals[i]) for i in range(nsph)]

    x = range(len(fitparams))
    labels = [fitparams[i] for i in range(len(fitparams))]
    oldparams = [oldparamdict[str(l)][fitparams[i]] for i in range(len(fitparams))]
    olderr = [oldparamdict[str(l)][fiterrparams[i]] for i in range(len(fiterrparams))]
    newparams = [newparamdict[str(l)][fitparams[i]] for i in range(len(fitparams))]
    newerr = [newparamdict[str(l)][fiterrparams[i]] for i in range(len(fiterrparams))]
    residual = (np.array(newparams)-np.array(oldparams))
    reserr = [np.sqrt(olderr[i]**2+newerr[i]**2) for i in range(len(olderr))]
    percent = (np.array(newparams)-np.array(oldparams))/np.array(oldparams)
    # percenterr = [100*np.sqrt((newparams[i]*olderr[i])**2/(oldparams[i])**4+(newerr[i]**2)/(oldparams[i]**2)) for i in range(len(olderr))]
    percenterr = [abs(percent[i])*np.sqrt(((newerr[i])/(newparams[i]))**2+((olderr[i])/(oldparams[i]))**2) for i in range(len(olderr))]

    fig, axarr = plt.subplots(2, sharex=True, figsize=(17,10))
    tPars = {'fontsize':16}
    plt.setp(axarr, xticks=x, xticklabels=labels)
    # axarr[0].errorbar(x,residual,yerr=reserr,linestyle="None", marker=".", markersize=10)
    axarr[0].hlines(0,-1,len(x),linestyle='-.')
    axarr[0].errorbar(x,newparams,yerr=newerr,linestyle="None", marker=".", markersize=13)
    #axarr[0].xticks(x,labels,**tPars)
    axarr[0].set_xlim(-1,len(fitparams))
    axarr[0].set_title(r'Old vs. new SH coefficients for $l_{}={}$'.format('{max}',l),**tPars)
    axarr[0].set_ylabel(r'SciPy Method Amplitudes',**tPars)
    # axarr[0].set_ylabel(r'Difference in Amplitudes',**tPars)
    axarr[0].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    axarr[1].errorbar(x,percent,yerr=percenterr,linestyle="None", marker=".", markersize=13)
    axarr[1].set_ylabel(r'Percent Difference in Amplitudes',**tPars)
    axarr[1].set_ylim(-0.1,0.1)
    plt.savefig('/home/jbourbeau/public_html/figures/2Dfit/SHcomparison_lmax{}_RIchi2.png'.format(l), dpi=300, bbox_inches='tight')



    # plt.show()
