#!/usr/bin/env python

import subprocess

from dataFunctions import getConfigs

if __name__ == "__main__":

    exList = []

    ##=======================================================================
    ## IceTop options

    configs = getConfigs('IT')
    filters = ['NotSTA8','STA8']
    methods = ['sid','anti','solar','ext']

    # All configs : 24H : NotSTA8|STA8 : sid|anti|solar|ext
    params  = [[c,f,m] for c in configs for f in filters for m in methods]
    for c, f, m in params:
        exList += ['./maker.py -c %s -f %s -m %s -n 20' % (c,f,m)]

    # All configs : 04H : NotSTA8|STA8 : sid
    params = [[c,f] for c in configs for f in filters]
    for c, f in params:
        exList += ['./maker.py -c %s -f %s -t 4 -n 20' % (c,f)]

    ##=======================================================================
    ## IceCube options

    configs = getConfigs('IC')
    methods = ['sid','anti','ext','solar']

    # All configs : 24H : sid|anti|solar|ext
    params = [[c, m] for c in configs for m in methods]
    for c, m in params:
        exList += ['./maker.py -c %s -m %s' % (c,m)]

    # All configs : 24H : sid : energy bins
    for c in configs:
        exList += ['./maker.py -c %s --ebins' % c]

    # All configs : 04H : sid
    for c in configs:
        exList += ['./maker.py -c %s -t 4' % c]

    for ex in exList:
        ex = ex.split(' ')
        subprocess.call(ex)
