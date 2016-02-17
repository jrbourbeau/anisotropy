#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import gammaln
import argparse

import load_sim
from powerFit import powerFit
from fakeSpec import fakeSpec
from llhtools import getEbins

# Compare shape of flux for each bin to ensemble average
def bayes(x1, x2):

    # Calculate totals
    n1 = float(np.sum(x1))
    n2 = float(np.sum(x2))
    output = (sum(gammaln(x1+1) + gammaln(x2+1) - gammaln(x1+x2+2))
            + gammaln(n1+n2+2) - gammaln(n1+1) - gammaln(n2+1))

    return output

def bayesTable(s, emin=6.2):

    t = np.log10(s['MC_energy'])
    r = np.log10(s['ML_energy'])
    pf_ratio = s['fLLH'] - s['pLLH']
    colorDict = {'P':'b','Fe':'r'}

    ebins = np.arange(6.25, 9.501, 0.25)
    llhbins = np.linspace(-15,15,100)

    for minE, maxE in zip(ebins[:-1], ebins[1:]):
        c0 = (t >= minE) * (t < maxE)
        c0 *= (r >= emin)
        c0 *= s['cuts']['llh']
        h = {}
        nevents = 0
        for true_comp in ['P','Fe']:
            c1 = c0 * (s['comp'] == true_comp)
            if c1.sum()==0:
                continue
            h[true_comp] = np.histogram(pf_ratio[c1], bins=llhbins)[0]
            nevents += c1.sum()

        print '%s %s %s %i' % (minE, maxE, nevents, bayes(h['P'], h['Fe']))


""" Goal: add iron events to proton sample until bayes factor reaches 100(?)"""
def stefan(s, cut=False, emin=6.2, out=False, spec=None):

    # Starting parameters
    r = np.log10(s['ML_energy'])
    c0 = (r >= emin)
    if cut:
        c0 *= s['cuts']['llh']
    llhbins = np.linspace(-10,10,400)
    # Create power spectrum cut
    if spec != None:
        temp = {'P':powerFit(spec), 'Fe':powerFit(spec)}
        specCut = fakeSpec(s, c0, temp, emin=emin)
        c0 *= specCut

    # Apply emin cut early
    pcut = (s['comp'] == 'P')[c0]
    fcut = (s['comp'] == 'Fe')[c0]
    pf_ratio = (s['fLLH'] - s['pLLH'])[c0]

    # Start with a pure proton distribution
    pdist = np.histogram(pf_ratio[pcut], bins=llhbins)[0]
    nproton = pdist.sum()
    print nproton

    niter = 100
    yvals = np.zeros((50,niter))
    for i in range(50):
        n = float(i) / 100
        niron = n * nproton / (1 - n)
        for n in range(100):
            fvals = np.random.choice(pf_ratio[fcut], niron)
            fdist = np.histogram(fvals, bins=llhbins)[0]
            mydist = pdist + fdist
            yvals[i][n] = bayes(pdist, mydist)

    x = range(50)
    y = np.median(yvals, axis=1)
    yup = np.percentile(yvals, 84, axis=1) - y
    ydn = y - np.percentile(yvals, 16, axis=1)
    yx = [np.log(100) for i in x]

    fig, ax = plt.subplots()
    ax.errorbar(x, y, yerr=[yup,ydn], fmt='.')
    ax.plot(x, yx)
    ax.set_xlabel(r'Iron in Sample [$\%$]', fontsize=14)
    ax.set_ylabel(r'Bayes Factor [$\ln(B_{21})$]', fontsize=14)
    #print yvals[0], yvals[-1]
    if out:
        plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.show()


""" Look for separation in true composition using the llh ratio """
def llhratio(s, emin=6.2, out=False):

    t = np.log10(s['MC_energy'])
    r = np.log10(s['ML_energy'])
    pf_ratio = s['fLLH'] - s['pLLH']
    colorDict = {'P':'b','Fe':'r'}

    llhbins = np.linspace(-5,5,100)

    fig, ax = plt.subplots()
    #ax.set_title('Likelihood Ratio', fontsize=18)
    ax.set_xlabel(r'Log-Likelihood Difference ($f_{\mathrm{LLH}}-p_{\mathrm{LLH}}$)', fontsize=16)
    ax.set_ylabel('Counts', fontsize=16)

    minE, maxE = 6.5, 6.75
    c0 = (t >= minE) * (t < maxE)
    c0 *= (r >= emin)
    c0 *= s['cuts']['llh']
    nevents = float(c0.sum())
    weights = np.array([1/nevents for i in c0])
    for true_comp in ['P','Fe']:
        c1 = c0 * (s['comp'] == true_comp)
        if c1.sum()==0:
            continue
        ax.hist(pf_ratio[c1], bins=llhbins, histtype='step',
                label=true_comp, linewidth=2, color=colorDict[true_comp])
        #        weights=weights[c1])
        # Custom text
        txtLabel = '%.2f-%.2f' % (minE, maxE)
        ax.set_xlim(-6,6)
        #ax.set_ylim(0,0.025)

    if out:
        plt.savefig(out)
    plt.show()


""" Look for separation in true composition using the llh ratio """
def llhratio_orig(s, n=4, emin=4, out=False):

    t = np.log10(s['MC_energy'])
    r = np.log10(s['ML_energy'])
    pf_ratio = s['fLLH'] - s['pLLH']
    colorDict = {'P':'b','Fe':'r'}

    ebins = getEbins(reco=True)
    ebins = ebins[ebins >= emin][::n]
    llhbins = np.linspace(-5,5,100)

    nrow, ncol = 3, 5
    #fig = plt.figure(figsize=(17,17))
    fig, axs = plt.subplots(nrow, ncol, sharex=True, sharey=True,
            figsize=(10,10))
    #ax = fig.add_subplot(111)
    #ax.set_title('Likelihood Ratio', fontsize=18)
    #ax.set_xlabel('Log-Likelihood Difference (fLLH-pLLH)', fontsize=16)
    #ax.set_ylabel('Counts', fontsize=16)

    for i, (minE, maxE) in enumerate(zip(ebins[:-1], ebins[1:])[:-1]):
        ri, ci = i/ncol, i%ncol
        c0 = (t >= minE) * (t < maxE)
        c0 *= (r >= emin)
        c0 *= s['cuts']['llh']
        nevents = float(c0.sum())
        weights = np.array([1/nevents for i in c0])
        for true_comp in ['P','Fe']:
            c1 = c0 * (s['comp'] == true_comp)
            if c1.sum()==0:
                continue
            axs[ri,ci].hist(pf_ratio[c1], bins=llhbins, histtype='step',
                    label=true_comp, linewidth=2, color=colorDict[true_comp],
                    weights=weights[c1])
            # Custom text
            txtLabel = '%.2f-%.2f' % (minE, maxE)
            axs[ri,ci].text(0.5, 0.85, txtLabel, ha='center',
                    transform=axs[ri,ci].transAxes)
            axs[ri,ci].set_xlim(-6,6)
            axs[ri,ci].set_ylim(0,0.025)

    # Text for axes
    fig.text(0, 0.5, '% of events in bin', va='center', rotation='vertical')
    fig.text(0.5, 0, 'Difference in likelihood (fLLH - pLLH)', ha='center')
    #ax.legend()

    if out:
        plt.savefig(out)
    plt.show()


if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument('-o', '--out', dest='out')
    args = p.parse_args()

    s = load_sim.load_sim()
    stefan(s, cut=True, emin=6.2, out=args.out, spec=-2.7)


