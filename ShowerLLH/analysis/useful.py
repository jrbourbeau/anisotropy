#!/usr/bin/env python

from numpy import *
from matplotlib.path import Path
import glob


##=========================================================================##
## Energy binning

""" Extract energy bins """
def getEbins():
    n_prefix = '/net/user/fmcnally/ShowerLLH/'
    binDict = load(n_prefix + 'resources/ShowerLLH_bins.npy')
    binDict = binDict.item()
    return binDict['Ebins'][44:-30]

""" Get average value of each energy bin """
def getEmids():
    Ebins = getEbins()
    Emids = log10((10**Ebins[1:] + 10**Ebins[:-1])/2.)
    return Emids


##=========================================================================##
## Cuts for use on data                          

""" Vertex finder for use with geometric cut """
def getVertex(d, config='IT73', tank_x=None, tank_y=None):

    # Define the sides by the tanks they contain, array built bottom up, L->R
    # Note: there appears to only be one tank showing at station 39
    side = {}
    if config == 'IT73':
        side[0] = array([0,1,2,3,4,5,6,7,8,9])
        side[1] = array([8,9,20,21,34,35,50,51,67,68,87,88])
        side[2] = array([87,88,105,106,121,122,135,136])
        side[3] = array([131,132,133,134,135,136])
        side[4] = array([131,132,143,144])
        side[5] = array([137,138,139,140,141,142,143,144])
        side[6] = array([69,70,89,90,107,108,123,124,137,138])
        side[7] = array([0,1,10,11,22,23,36,37,52,53,69,70])
    if config == 'IT81':
        side[0] = array([0,1,2,3,4,5,6,7,8,9,10,11])
        side[1] = array([10,11,24,25,40,41,58,59,77,78,97,98])
        side[2] = array([97,98,115,116,131,132,145,146])
        side[3] = array([141,142,143,144,145,146])
        side[4] = array([141,142,153,154])
        side[5] = array([147,148,149,150,151,152,153,154])
        side[6] = array([60,61,79,80,99,100,117,118,133,134,147,148])
        side[7] = array([0,1,12,13,26,27,42,43,60,61])
    '''
    if not tank_x or not tank_y:
        resourcedir = '/net/user/fmcnally/ShowerLLH/resources/'
        xy = load(resourcedir + config+'_tankpos.npy')
        tank_x = xy[0]
        tank_y = xy[1]
    '''
    # Find slopes and y-intercepts
    nside = len(side)
    m, b = zeros((2,nside))
    for i in range(nside):
        m[i], b[i] = polyfit(tank_x[side[i]], tank_y[side[i]], 1)

    # Intersection function
    def x_intersect(m1, b1, m2, b2):
        x = (b2-b1)/(m1-m2)
        return x

    # d = distance moved perpendicular from original line position
    theta = arctan(m)
    b1 = b + sign(b)*d/cos(theta)
    x = array([x_intersect(m[i-1], b1[i-1], m[i], b1[i]) for i in range(nside)])
    y = m*x + b1

    return x, y


""" Geometric containment cut """
def inPoly(xx, yy, d, config='IT73'):

    # Create shape of detector
    vertx, verty = getVertex(d, config=config)
    if vertx[-1]!=vertx[0]:
        vertx = append(vertx, vertx[0])
        verty = append(verty, verty[0])
    xyverts  = array([vertx, verty]).transpose()
    # Corresponding codes
    codes = [Path.LINETO for i in xyverts]
    codes[0]  = Path.MOVETO
    codes[-1] = Path.CLOSEPOLY
    path = Path(xyverts, codes)

    xypoints = array([xx,yy]).transpose()
    tf = path.contains_points(xypoints)

    return tf



##=========================================================================##
## Other useful functions


""" Return an array of the number of counts per energy bin """
def Nfinder(a, cut, w=False):
    Ebins = getEbins()
    if not w:
        Npassed, bins = histogram(a[cut], Ebins)
    else:
        Npassed = zeros(len(Ebins)-1)
        for wkey in w.keys():
            weight = float(wkey[1:])
            temp_n, bins = histogram(a[cut*w[wkey]], Ebins)
            Npassed += weight * temp_n
    return Npassed.astype('float')


""" Return an array of the error per energy bin """
def Efinder(a, cut, w=False):
    Ebins = getEbins()
    if not w:
        x = a[cut]
        n, bins = histogram(x, Ebins)
    else:
        n = zeros(len(Ebins)-1)
        for wkey in w.keys():
            weight = float(wkey[1:])
            temp_n, bins = histogram(a[cut*w[wkey]], Ebins)
            n += weight**2 * temp_n
    return sqrt(n)



#############################################################################
# Less useful and/or just used for plotting tests
#############################################################################


""" Simple smoothing function for unfolded spectrum """
def smoother(a):

    temp = zeros(len(a))
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

