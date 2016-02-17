#!/usr/bin/env python

##============================================================================
## Temporary hack to limit the energy bins and likelihood tables to 5-9.5
## Run after TableConverter
##============================================================================

import myGlobals as my
import numpy as np

if __name__ == "__main__":

    my.setupShowerLLH(verbose=False)
    f = '%s/LLHTables_logdist.npy' % my.llh_resource
    d = np.load(f)
    d = d.item()

    d['bins'][0][1] = np.linspace(5, 9.5, 91)
    for key in d['llhtables'].keys():
        d['llhtables'][key] = d['llhtables'][key][20:]

    np.save(f, d)
