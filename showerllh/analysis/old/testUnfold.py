#!/usr/bin/env python

from numpy import *
import sim, unfold, sys, powerFit
sys.path.append('/home/fmcnally/ShowerLLH/useful/')
from LLHFunctions import *
import matplotlib.pyplot as plt

##============================================================================
## Useful functions

def crapE(true):
    mu = 10**true
    sigma = 0.5 * 10**true
    reco = log10(random.normal(mu, sigma))
    return reco


def crapComp(true):
    sample = random.uniform(size=len(true))
    fakep = (sample < 0.5)
    fakef = (sample > 0.5)
    return fakep, fakef


def crapProb(s, cut):

    q = {}
    true = log10(s['MC_energy'])
    reco = crapE(true)
    # 50% smearing gaussian can cause energies to have negative values
    nancut = (reco==reco)
    comp = {}
    comp['p'], comp['f'] = crapComp(true)

    # Calculate probabilities
    erTable = [[e, r] for e in ['F','P'] for r in ['f','p']]
    for e, r in erTable:
        c0 = cut * s[e] * comp[r] * nancut
        RvT, x, y = histogram2d(reco[c0], true[c0], bins=(Ebins,Ebins))
        q['RvT_'+e+r] = asarray(RvT, dtype='float')

    # Normalize
    Ntot = {}
    Ntot['F'] = sum(q['RvT_Ff'] + q['RvT_Fp'], axis=0)
    Ntot['P'] = sum(q['RvT_Pf'] + q['RvT_Pp'], axis=0)
    for e, r in erTable:
        q['RvT_'+e+r] /= Ntot[e]

    # Store in probability dictionary
    # p['Rr|Tt'][i][j] = probability comp t in bin j is reconstructed as r in i
    #  - sum each column to get 1
    prob = {'Rf|Tf':q['RvT_Ff'], 'Rf|Tp':q['RvT_Pf']}
    prob.update({'Rp|Tp':q['RvT_Pp'], 'Rp|Tf':q['RvT_Fp']})

    # Setup so N_true for lowest energy is >= N_true for highest energy
    N_max = Nfinder(true, cut)[-1]
    min_bin = nonzero(Nfinder(true, cut) > N_max)[0][0]
    for key in prob.keys():
        prob[key] = prob[key][:,min_bin:]

    return prob


##===========================================================================##
## Applying unfolding to simulation

## Plot the counts after unfolding ##
def counts(s, niter, sDict={'':[-0.25, [[7.5, -0.75]]]},
                        smooth=True, spl=False, diff=False):

    t = log10(s['MC_energy'])
    r = crapE(t)
    cutName = 'llh'
    nancut = (r==r)
    cut = s['cuts'][cutName] * nancut

    f = plt.figure()
    ax1 = f.add_subplot(1,1,1)
    ax1.set_title('Energy Spectum using '+cutName+' cut ('+str(niter)+' iterations)')
    ax1.set_xlabel('Log10(Energy/GeV)')
    ax1.set_ylabel('Counts')

    # Load probability tables
    print 'Getting probability tables...'
    p = crapProb(s, cut)
    st = len(sim.Emids) - len(p['Rf|Tf'][0])
    Emids = sim.Emids[st:]

    # Option for splining
    if spl:
        for key in p.keys():
            p[key] = 10**(spline.spline(s, p[key], nk=2, npoly=3))
        print sum(p['Rf|Tp']+p['Rp|Tp'], axis=0)
        print sum(p['Rf|Tf']+p['Rp|Tf'], axis=0)

    # Create our toy MC spectrum
    temp = {}
    for key in sDict.keys():
        s0 = sDict[key][0]
        try:
            sTable = sDict[key][1]
        except IndexError:
            sTable=False
        temp[key] = powerFit.powerFit(s0, sTable=sTable, st=st)
    specCut = sim.fakeSpec(s, cut, temp)
    #specCut = array([True for i in range(len(specCut))])

    # Create starting arrays
    N_passed = Nfinder(r, cut*specCut).sum()
    comp = {}
    comp['p'], comp['f'] = crapComp(t)
    N_proton = Nfinder(r, cut*comp['p']*specCut)
    N_iron   = Nfinder(r, cut*comp['f']*specCut)

    # Get starting probabilities
    p['Rp'] = N_proton / N_passed
    p['Rf'] = N_iron / N_passed
    p['Tf'], p['Tp'] = [powerFit.powerFit(-2.7, st=st) for i in range(2)]

    # Get relative errors due to unfolding
    # Due to efficiency
    effarea, relerr = sim.getEff(s, cut, smooth=smooth)
    relerr = relerr[st:]

    # Bayesian unfolding
    for i in range(niter):
        p['Tf'], p['Tp'] = unfold.unfold(p)
        # Smooth prior before next iteration (except last time)
        if i < niter-1:
            p['Tf'] = smoother(p['Tf'])
            p['Tp'] = smoother(p['Tp'])

    # Find bin values and errors
    F_Nunfold = p['Tf'] * N_passed
    P_Nunfold = p['Tp'] * N_passed
    All_Nunfold = F_Nunfold + P_Nunfold
    ## NOTE: I don't think you can use the relative errors like this ##
    F_relerr = sqrt(1/F_Nunfold + relerr**2)
    P_relerr = sqrt(1/P_Nunfold + relerr**2)
    F_err = F_Nunfold * F_relerr
    P_err = P_Nunfold * P_relerr
    All_err = sqrt(F_err**2 + P_err**2)

    #ax1.errorbar(Emids, F_Nunfold, yerr=F_err, fmt='r.', label='Iron')
    #ax1.errorbar(Emids, P_Nunfold, yerr=P_err, fmt='b.', label='Proton')

    # Plot true spectrum
    MC_N = Nfinder(log10(s['MC_energy']), cut*specCut)[st:]
    MC_F = Nfinder(log10(s['MC_energy']), cut*specCut*s['F'])[st:]
    MC_P = Nfinder(log10(s['MC_energy']), cut*specCut*s['P'])[st:]
    #ax1.plot(Emids, MC_F, 'rx', label='MC_F')
    #ax1.plot(Emids, MC_P, 'bx', label='MC_P')

    # plot original (not unfolded) spectrum
    O_N = (N_proton + N_iron)[st:]
    O_relerr = sqrt(1/O_N + relerr**2)
    O_err = O_N * O_relerr

    if not diff:
        ax1.errorbar(Emids, O_N, yerr=O_err, fmt='k.', label='Original')
        ax1.plot(Emids, MC_N, 'rx', label='MC')
        ax1.errorbar(Emids, All_Nunfold, yerr=All_err, fmt='g.', label='Unfold')
        ax1.set_yscale('log')

    if diff:
        ax1.plot(Emids, (O_N - MC_N), 'k', label='Original - MC')
        ax1.plot(Emids, (All_Nunfold - MC_N), 'r', label='Unfold - MC')

    ax1.legend(loc='lower left')
    #ax1.set_ylim((10**(-1), 10**(4)))
    plt.show()


