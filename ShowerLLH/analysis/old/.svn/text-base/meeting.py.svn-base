#!/usr/bin/env python

from numpy import *
import sys, pickle, sim, unfold, powerFit, subprocess
sys.path.append('/home/fmcnally/ShowerLLH/useful')
from LLHFunctions import *
from data import NumToFlux as n2f
from data import smoother
from load_data import load_data

def histWriter(d, s, mp='All'):

    ##=========================================================================
    ## Basic setup

    # Starting information
    r = log10(d['ML_energy'])
    wList = [d['w1'], d['w24'], d['w8']]
    cutName = 'llh'
    dcut = d['cuts'][cutName]
    scut = s['cuts'][cutName]

    monthList = [mp]
    if mp == 'All':
        monthList = ['201006', '201007', '201008', '201012', '201101']

    # Get probabilities for unfolding
    fl = open('probabilities.pkl', 'rb')
    p = pickle.load(fl)
    fl.close()
    st = len(sim.Emids) - len(p['Rf|Tf'][0])
    Emids = sim.Emids[st:]

    # Relative errors
    fl = open('unfold_err.pkl', 'rb')
    chi2, unrel = pickle.load(fl)
    fl.close()
    effarea, effrel = sim.getEff(s, scut, smooth=True)
    effarea, effrel = effarea[st:], effrel[st:]


    ##=========================================================================
    ## Make histograms

    q = {}
    q['Emids'] = Emids
    q['flux'], q['err'], q['chi2'] = {},{},{}
    d['cuts']['none'] = array([True for i in range(len(dcut))])
    nameList = ['none', 'left', 'right', 'odds', 'evens']
    cutList = [d['cuts'][name] for name in nameList]
    niter = 30

    for k in range(len(cutList)):

        # Get counts
        c0 = dcut * cutList[k]
        N_passed = Nfinder(r, c0, wList).sum()
        N_proton = Nfinder(r, c0*d['p'], wList)
        N_iron   = Nfinder(r, c0*d['f'], wList)
        name = nameList[k]

        # Get starting probabilities
        p['Rp'] = N_proton / N_passed
        p['Rf'] = N_iron / N_passed
        p['Tp'], p['Tf'] = [powerFit.powerFit(-2.7, st=st) for i in range(2)]
        q['chi2'][name] = []

        # Bayesian unfolding
        Tf, Tp = p['Tf'], p['Tp']
        for i in range(niter):

            # Calculate unfolded fluxes, errors, and chi sqaured values
            p['Tf'], p['Tp'] = unfold.unfold(p)
            old = Tf + Tp
            new = p['Tf'] + p['Tp']
            q['chi2'][name].append( 1./(len(Emids)+1) * 
                    sum(((new-old) / unrel[i])**2))

            counts, flux, relerr, err = {},{},{},{}
            counts['A'] = (p['Tf']+p['Tp']) * N_passed
            counts['F'] = p['Tf'] * N_passed
            counts['P'] = p['Tp'] * N_passed
            for key in counts.keys():
                relerr[key] = sqrt(1/counts[key] + effrel**2 + unrel[i]**2)
                flux[key] = n2f(counts[key], effarea, monthList, st)
                err[key] = flux[key] * relerr[key]
                # Write to dictionary
                q['flux'][name+'_'+key+'_'+str(i+1)] = flux[key]
                q['err'][name+'_'+key+'_'+str(i+1)] = err[key]

            if i < niter-1:
                p['Tf'] = smoother(p['Tf'])
                p['Tp'] = smoother(p['Tp'])

        # Original (not unfolded) spectrum
        O_counts = (N_proton + N_iron)[st:]
        O_relerr = sqrt(1/O_counts + effrel**2)
        q['flux'][name+'_O'] = n2f(O_counts, effarea, monthList, st)
        q['err'][name+'_O'] = q['flux'][name+'_O'] * O_relerr


    ##=========================================================================
    ## Write to file

    print 'Writing to file...'
    fl = open('collab/'+mp+'_hists.pkl', 'wb')
    pickle.dump(q, fl)
    fl.close()


if __name__ == "__main__":

    # Load simulation
    print 'Loading simulation...'
    vars = {}
    execfile('load_sim.py', vars)
    s = vars['s']

    #for month in ['201006', '201007', '201008', '201012', '201101', 'All']:
    for month in ['All']:
        # Load data
        print 'Loading', month
        d, trash = load_data(month)
        histWriter(d, s, mp=month)








