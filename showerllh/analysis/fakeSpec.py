#!/usr/bin/env python

##############################################################################
# Tools for creating a fake spectrum from simulation
##############################################################################

import numpy as np
from llhtools import *
from useful import getMids

def fakeSpec(s, cut, testSpec, reco=True, emin=None):

    print 'Creating fake spectrum...'
    t = np.log10(s['MC_energy'])
    Ebins = getEbins(reco=reco)
    Emids = getMids(Ebins)
    if emin != None:
        Emids = Emids[Emids > emin]
        Ebins = Ebins[-len(Emids)-1:]
        for e in testSpec.keys():
            testSpec[e] = testSpec[e][-len(Emids):]

    #st = len(Emids) - len(testSpec[testSpec.keys()[0]])
    samples = np.zeros(len(t), dtype='bool')
    passindex = []
    #binList = np.digitize(t, Ebins[st:]) - 1
    binList = np.digitize(t, Ebins) - 1

    # testSpec should contain spectrum shape for proton and iron
    for e in testSpec.keys():

        compCut = np.array([True for i in range(len(cut))])
        if e in ['P','He','O','Fe']:
            compCut = s['comp'] == e

        c0 = cut * compCut
        #N_passed = Nfinder(t, c0)[st:]
        N_passed = np.histogram(t[c0], bins=Ebins)[0]

        # Fit the distribution to our available data
        diff = N_passed / testSpec[e]
        shift = N_passed[diff.argmin()] / testSpec[e][diff.argmin()]
        testSpec[e] *= shift
        # Wiggle within poisson errors
        testSpec[e] = np.around(testSpec[e])
        for i in range(len(testSpec[e])):
            if i != diff.argmin():
                testSpec[e][i] = np.random.poisson(testSpec[e][i])

        # From distribution, pull out events to create spectrum
        #for i in range(len(Emids[st:])):
        for i in range(len(Emids)):
            # Get events in desired energy bin
            ecut = (binList == i)
            mc = (ecut * c0).sum()
            # Build array of randomly selected indexes
            temp = np.zeros(mc, dtype='bool')
            temp[:testSpec[e][i]] = True
            np.random.shuffle(temp)
            # Find indexes of events that passed
            passindex.append(np.where(ecut*c0 == True) * temp)

    # Fill samples with events that passed
    for index in passindex:
        samples[index] = True

    return samples

