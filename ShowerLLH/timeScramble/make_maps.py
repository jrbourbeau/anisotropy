#!/usr/bin/env python

from numpy import *
import sys, time
from timeScramble import timeScramble
sys.path.append('/home/zgriffith/ShowerLLH/analysis')
import load_data, llh_cut


def maker(d, cut, nInt, method='equatorial', out=False):

    t0 = time.time()
    data = {}
    # Apply cuts
    for key in ['mjd', 'ShowerPlane_zenith', 'ShowerPlane_azimuth']:
        data[key] = d[key][cut]
    for key in d['weights'].keys():
        data[key] = d['weights'][key][cut]

    timeScramble(data, nInt, method=method, out=out)
    print 'Scrambling time:', time.time()-t0


def ecuts_maker(d, nInt, method='equatorial', out=False):

    # Setup
    llh_energy = log10(d['ML_energy'])
    c0 = d['cuts']['llh'] * d['weights']['w1']
    #Ebins = arange(5, 9.75, 0.5)
    Ebins = range(5,10)

    # Loop through energy bins
    for i in range(len(Ebins)-1):
        e_cut = (llh_energy>=Ebins[i]) * (llh_energy<Ebins[i+1])
        cut = c0 * e_cut
        if cut.sum() == 0:
            continue
        print 'Working on %s to %s' % (Ebins[i], Ebins[i+1])
        if out:
            temp_out = '%s_%sGeV_%sGeV' % (out, Ebins[i], Ebins[i+1])
        maker(d, cut, nInt, method=method, out=temp_out)


def nstation_maker(d, nInt, method='equatorial', out=''):

    # Setup - bins roughly equal counts from data
    binList = array([0, 4, 5, 6, 7, 8, 10, 15, 99])
    c0 = d['cuts']['llh']

    for i in range(len(binList) - 1):
        n_min, n_max = binList[i], binList[i+1]
        n_cut = (d['NStations'] >= n_min) * (d['NStations'] < n_max)
        cut = c0 * n_cut
        temp_out = '%s_%02ista_%02ista' % (out, n_min, n_max)
        maker(d, cut, nInt, method=method, out=temp_out)


def llh_maker(d, nInt, method='equatorial', out=''):

    # Setup
    c0 = d['cuts']['llh']
    for comp in ['p','f']:
        #cut = c0 * d['cuts'][comp+'LLH']
        cut = c0 * d[comp]
        temp_out = '%s_%sLLH' % (out, comp)
        maker(d, cut, nInt, method=method, out=temp_out)


if __name__ == "__main__":

    config = sys.argv[1]
    ymList = sys.argv[2:]
    nInt = 24
    method = 'equatorial'
    mapDir = '/net/user/zgriffith/ShowerLLH/maps/raw'

    for yyyymm in ymList:

        d = load_data.load_data(config, yyyymm)
        out = '%s/%s_%s_%sH' % (mapDir, config, yyyymm, nInt)
        #cut = d['cuts']['llh']

        #ecuts_maker(d, nInt, method=method, out=out)
        #nstation_maker(d, nInt, method=method, out=out)
        #llh_maker(d, nInt, method=method, out=out)
        #cut = d['cuts']['llh']
        #maker(d, cut, nInt, method=method, out=out)
        #cut = d['cuts']['llh'] * d['weights']['w1']
        #maker(d, cut, nInt, method=method, out=out+'_STA8')
        #cut = d['cuts']['llh'] * d['weights']['w1.5']
        #maker(d, cut, nInt, method=method, out=out+'_STA3ii')
        #cut = d['cuts']['llh'] * d['weights']['w3']
        #maker(d, cut, nInt, method=method, out=out+'_STA3')
        cut = d['cuts']['llh'] * logical_not(d['weights']['w3'])
        maker(d, cut, nInt, method=method, out=out+'_NotSTA8')


