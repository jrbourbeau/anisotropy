#!/usr/bin/env python

###############################################################################
## Tests the array mathematics, proving that our method is identical to the
## slower method of looping over a for loop
###############################################################################

from numpy import *
from random import random


def testP(n, m):

    # Create normalized probability tables for n true bins, m reco bins
    rptp = array([array([random() for i in range(n)]) for j in range(m)])
    rptf = array([array([random() for i in range(n)]) for j in range(m)])
    rftp = array([array([random() for i in range(n)]) for j in range(m)])
    rftf = array([array([random() for i in range(n)]) for j in range(m)])
    ptot = sum(rptp+rftp, axis=0)
    ftot = sum(rptf+rftf, axis=0)
    rptp /= ptot
    rftp /= ptot
    rptf /= ftot
    rftf /= ftot

    # Make sample starter probability
    tp = array([random() for i in range(n)])
    tf = array([random() for i in range(n)])
    ptot = sum(tp+tf)
    tp /= ptot
    tf /= ptot

    # Calculate probability matrices the indexed (slow) way
    tprp, tprf, tfrp, tfrf = zeros((4, n, m))
    for i in range(n):
        for j in range(m):
            tprp[i][j] = (rptp[j][i]*tp[i] / 
                sum([rptp[j][k]*tp[k] + rptf[j][k]*tf[k] for k in range(n)]))
            tprf[i][j] = (rftp[j][i]*tp[i] / 
                sum([rftp[j][k]*tp[k] + rftf[j][k]*tf[k] for k in range(n)]))
            tfrp[i][j] = (rptf[j][i]*tf[i] / 
                sum([rptp[j][k]*tp[k] + rptf[j][k]*tf[k] for k in range(n)]))
            tfrf[i][j] = (rftf[j][i]*tf[i] / 
                sum([rftp[j][k]*tp[k] + rftf[j][k]*tf[k] for k in range(n)]))

    # Calculate probability matrices our way
    p = {}
    p['Tp|Rp'] = (rptp*tp).transpose() / sum(rptf*tf + rptp*tp, axis=1)
    p['Tp|Rf'] = (rftp*tp).transpose() / sum(rftf*tf + rftp*tp, axis=1)
    p['Tf|Rp'] = (rptf*tf).transpose() / sum(rptf*tf + rptp*tp, axis=1)
    p['Tf|Rf'] = (rftf*tf).transpose() / sum(rftf*tf + rftp*tp, axis=1)

    print 'tprp =', tprp
    print 'ours =', p['Tp|Rp']
    print 'tprf =', tprf
    print 'ours =', p['Tp|Rf']
    print 'tfrp =', tfrp
    print 'ours =', p['Tf|Rp']
    print 'tfrf =', tfrf
    print 'ours =', p['Tf|Rf']


if __name__ == "__main__":

    probTest(3, 4)

