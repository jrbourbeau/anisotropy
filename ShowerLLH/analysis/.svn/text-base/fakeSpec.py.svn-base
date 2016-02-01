#!/usr/bin/env python

##############################################################################
# Tools for creating a fake spectrum from simulation
##############################################################################

from numpy import *
from useful import *

def fakeSpec(s, cut, testSpec):

    print 'Creating fake spectrum...'
    t = log10(s['MC_energy'])
    Ebins, Emids = getEbins(), getEmids()
    st = len(Emids) - len(testSpec[testSpec.keys()[0]])
    samples = zeros(len(t), dtype='bool')
    passindex = []
    binList = digitize(t, Ebins[st:]) - 1

    # testSpec should contain spectrum shape for proton and iron
    for e in testSpec.keys():

        compCut = array([True for i in range(len(cut))])
        if e in ['P','He','O','Fe']:
            compCut = s[e]

        c0 = cut * compCut
        N_passed = Nfinder(t, c0)[st:]

        # Fit the distribution to our available data
        diff = N_passed / testSpec[e]
        shift = N_passed[diff.argmin()] / testSpec[e][diff.argmin()]
        testSpec[e] *= shift
        # Wiggle within poisson errors
        testSpec[e] = around(testSpec[e])
        for i in range(len(testSpec[e])):
            if i != diff.argmin():
                testSpec[e][i] = random.poisson(testSpec[e][i])

        # From distribution, pull out events to create spectrum
        for i in range(len(Emids[st:])):
            # Get events in desired energy bin
            ecut = (binList == i)
            mc = (ecut * c0).sum()
            # Build array of randomly selected indexes
            temp = zeros(mc, dtype='bool')
            temp[:testSpec[e][i]] = True
            random.shuffle(temp)
            # Find indexes of events that passed
            passindex.append(where(ecut*c0 == True) * temp)

    # Fill samples with events that passed
    for index in passindex:
        samples[index] = True

    return samples

