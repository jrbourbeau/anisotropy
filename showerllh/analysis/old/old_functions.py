#!/usr/bin/env python

###############################################################################
# A collection of functions no longer used but potentially useful
###############################################################################


""" VEM to GeV conversion taken from Bakhtiyar's work with S125 """
def vem2gev2(s125, comp, zenith):
    a = cos(zenith)
    if comp == 'P':
        if a >= 0.95:
            p = [5.9977, 0.9589, -0.001]
        elif a >= 0.9:
            p = [6.0340, 0.9455, -0.0007]
        elif a >= 0.85:
            p = [6.0823, 0.9288, -0.0013]
        elif a >= 0.8:
            p = [6.1397, 0.9186, -0.0013]
        else:
            p = [0, 0, 0]
    elif comp == 'F':
        if a >= 0.95:
            p = [6.0893, 0.8695, 0.0143]
        elif a >= 0.9:
            p = [6.1525, 0.8516, 0.0164]
        elif a >= 0.85:
            p = [6.2241, 0.8427, 0.0147]
        elif a >= 0.8:
            p = [6.3066, 0.8372, 0.0134]
        else:
            p = [0, 0, 0]

    logE = p[2]*log10(s125)**2 + p[1]*log10(s125) + p[0]
    return logE

