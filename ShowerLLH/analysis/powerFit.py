#!/usr/bin/env python

from numpy import *
import matplotlib.pyplot as plt
import sim
from useful import getEmids
from scipy import optimize

def powerFit(s0, sTable=False, st=False):

    Emids = getEmids()
    if st:
        Emids = Emids[st:]
    Tf = (10**Emids)**s0
    if sTable:
        for e, spec in sTable:
            start = nonzero(Emids >= e)[0][0]
            shift = Tf[start] / ((10**Emids[start])**spec)
            Tf[start:] = shift * (10**Emids[start:])**spec

    sumprobs = Tf.sum() * 2
    Tf /= sumprobs

    return Tf


""" Fit n power laws to a given spectrum """
def pow_fit(x, y, sigmay, p0):

    # Start with the assumption of one break
    def powerLaw(x, s0, offset, e1, s1):
        testSpec = x**s0
        start = nonzero(x >= 10**e1)[0][0]
        shift = testSpec[start] / ((x[start])**s1)
        testSpec[start:] = shift * (x[start:])**s1
        testSpec *= offset
        return testSpec

    popt, pcov = optimize.curve_fit(powerLaw, x, y, p0, sigma=sigmay)
    print popt, pcov
    yfit = powerLaw(x, *popt)

    return yfit

def gaussFit(a, b, c, st=False):

    Emids = getEmids()
    if st:
        Emids = Emids[st:]
    Tf = a * exp(-(Emids-b)**2/(2*c**2))
    sumprobs = Tf.sum() * 2
    Tf /= sumprobs

    return Tf



if __name__ == "__main__":

    s0 = -2.7
    sTable = [[6, -3.0], [7.5, -2.7]]
    Tf, Tp = powerFit(s0, sTable, 0)

    Emids = getEmids()
    plt.plot(Emids, Tf, '.')
    plt.yscale('log')
    plt.show()

