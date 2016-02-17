#!/usr/bin/env python

import numpy as np
from matplotlib.path import Path
import glob

import myGlobals as my
my.setupShowerLLH(verbose=False)

##=========================================================================##
## Energy binning

""" Extract energy bins """
#def getEbins(bintype='logdist'):
#    binDict = np.load('%s/ShowerLLH_bins.npy' % my.llh_resource)
#    binDict = binDict.item()
#    bins = binDict['logdist']
#    binVals = [bins[i][0] for i in bins]
#    idx = binVals.index('E')
#    ebins = bins[idx][1]
#    return ebins
""" Temporary replacement """
def getEbins(reco=False):
    ebins = np.arange(4, 9.51, 0.05)
    if reco:
        ebins = ebins[20:]
    return ebins


""" Get average value of each energy bin """
def getEmids(bintype='logdist'):
    Ebins = getEbins(bintype=bintype)
    Emids = np.log10((10**Ebins[1:] + 10**Ebins[:-1])/2.)
    return Emids

##=========================================================================##
## Composition-specific information

def getComps(d, reco=True):
    if reco:
        recos = [key for key in d.keys() if 'ML_x' in key and len(key)==5]
        comps = [r[0] for r in recos]
        return comps
    else:
        trues = list(set(d['comp']))
        return trues

def getColor(e):
    colorDict = {'p':'b', 'h':'y', 'o':'c', 'f':'r', 'All':'k'}
    colorDict.update({'P':'b', 'He':'y', 'O':'c', 'Fe':'r'})
    return colorDict[e]



##=========================================================================##
## Cuts for use on data                          

""" Vertex finder for use with geometric cut """
def getVertex(d, config='IT73', tank_xy=None):

    ## IMPORTANT: inPoly will fail for IT59 - the last side moves the wrong way

    # Define the sides by the tanks they contain, array built bottom up, L->R
    # Note: there appears to only be one tank showing at station 39
    side = {}
    if config == 'IT59':
        side[0]  = np.array([0,1,2,3,4,5,6,7,8,9])
        side[1]  = np.array([8,9,18,19,28,29,38,39,47,48,61,62])
        side[2]  = np.array([61,62,77,78,93,94,107,108])
        side[3]  = np.array([107,108,105,106,103,104])
        side[4]  = np.array([103,104,115,116])
        side[5]  = np.array([115,116,113,114,111,112,109,110])
        side[6]  = np.array([109,110,95,96,79,80])
        side[7]  = np.array([79,80,63,64])
        side[8]  = np.array([63,64,65,66])
        side[9]  = np.array([65,66,49,50])
        side[10] = np.array([49,50,51,52])
        side[11] = np.array([51,52,40,41])
        side[12] = np.array([40,41,30,31,20,21,10,11,0,1])
    if config == 'IT73':
        side[0] = np.array([0,1,2,3,4,5,6,7,8,9])
        side[1] = np.array([8,9,20,21,34,35,50,51,67,68,87,88])
        side[2] = np.array([87,88,105,106,121,122,135,136])
        side[3] = np.array([131,132,133,134,135,136])
        side[4] = np.array([131,132,143,144])
        side[5] = np.array([137,138,139,140,141,142,143,144])
        side[6] = np.array([69,70,89,90,107,108,123,124,137,138])
        side[7] = np.array([0,1,10,11,22,23,36,37,52,53,69,70])
    if config == 'IT81':
        side[0] = np.array([0,1,2,3,4,5,6,7,8,9,10,11])
        side[1] = np.array([10,11,24,25,40,41,58,59,77,78,97,98])
        side[2] = np.array([97,98,115,116,131,132,145,146])
        side[3] = np.array([141,142,143,144,145,146])
        side[4] = np.array([141,142,153,154])
        side[5] = np.array([147,148,149,150,151,152,153,154])
        side[6] = np.array([60,61,79,80,99,100,117,118,133,134,147,148])
        side[7] = np.array([0,1,12,13,26,27,42,43,60,61])

    if tank_xy==None:
        import myGlobals as my
        my.setupShowerLLH(verbose=False)
        tank_xy = np.load('%s/tankpos.npy' % my.llh_resource)
        tank_xy = tank_xy.item()
        tank_xy = tank_xy[config]

    # Find slopes and y-intercepts
    nside = len(side)
    m, b = np.zeros((2,nside))
    tank_x, tank_y = np.transpose(tank_xy)
    for i in range(nside):
        m[i], b[i] = np.polyfit(tank_x[side[i]], tank_y[side[i]], 1)

    # Intersection function
    def x_intersect(m1, b1, m2, b2):
        x = (b2-b1)/(m1-m2)
        return x

    # d = distance moved perpendicular from original line position
    theta = np.arctan(m)
    b1 = b + np.sign(b)*d/np.cos(theta)
    x = np.array([x_intersect(m[i-1], b1[i-1], m[i], b1[i]) \
            for i in range(nside)])
    y = m*x + b1

    return x, y


""" Geometric containment cut """
def inPoly(xx, yy, d, config='IT73'):

    # Create shape of detector
    vertx, verty = getVertex(d, config=config)
    if vertx[-1]!=vertx[0]:
        vertx = np.append(vertx, vertx[0])
        verty = np.append(verty, verty[0])
    xyverts  = np.transpose([vertx, verty])
    # Corresponding codes
    codes = [Path.LINETO for i in xyverts]
    codes[0]  = Path.MOVETO
    codes[-1] = Path.CLOSEPOLY
    path = Path(xyverts, codes)

    xypoints = np.transpose([xx,yy])
    tf = path.contains_points(xypoints)

    return tf



##=========================================================================##
## Other useful functions


""" Return an array of the number of counts per energy bin """
def Nfinder(a, cut, w=None, bintype='logdist'):
    Ebins = getEbins()
    wtemp = None if w==None else w[cut]
    Npassed, bins = np.histogram(a[cut], Ebins, weights=wtemp)
    return Npassed.astype('float')


""" Return an array of the error per energy bin """
def Efinder(a, cut, w=None, bintype='logdist'):
    Ebins = getEbins()
    wtemp = None if w==None else w[cut]**2
    n, bins = np.histogram(a[cut], Ebins, weights=wtemp)
    return np.sqrt(n)



#############################################################################
# Less useful and/or just used for plotting tests
#############################################################################


""" Simple smoothing function for unfolded spectrum """
def smoother(a):

    temp = np.zeros(len(a))
    for i in range(len(a)):
        if i == 0:
            temp[i] = (2*a[i] + a[i+1]) / 3
        elif i == len(a)-1:
            temp[i] = (2*a[i] + a[i-1]) / 3
        else:
            temp[i] = (a[i-1] + 2*a[i] + a[i+1]) / 4
    return temp


""" Get the information for IT26 spectrum """
def getBakh():
    file = genfromtxt('it73.dat')
    logE = file[:,0]
    Nevents = file[:,1]
    flux = file[:,2]
    sigma = file[:,3]
    syserr_upper = file[:,4]
    syserr_lower = file[:,5]
    return [logE, Nevents, flux, sigma, syserr_upper, syserr_lower]

