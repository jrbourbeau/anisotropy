#!/usr/bin/env python

import numpy as np
import tables

from icecube import astro

if __name__ == "__main__":


    infile = '/data/user/fmcnally/showerllh/IT81-II_data/files/DataLLH_20120630_logdist.hdf5'

    t = tables.openFile(infile)
    t1 = astro.Time()

    print 'Loading information...'
    mjd_day = t.root.I3EventHeader.col('time_start_mjd_day')
    mjd_sec = t.root.I3EventHeader.col('time_start_mjd_sec')
    mjd_ns  = t.root.I3EventHeader.col('time_start_mjd_ns')

    mjd = np.zeros(len(mjd_day), dtype=np.float64)
    for i in range(len(mjd_day)):
        day = int(mjd_day[i])
        sec = int(mjd_sec[i])
        ns  = int(mjd_ns[i])
        t1.SetTime(day, sec, ns)
        mjd[i] = t1.GetMJD()

    t.close()

