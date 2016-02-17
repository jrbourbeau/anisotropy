#!/usr/bin/env python

""" Calculates time difference between 1st and last event"""
def duration(d):

    #from icecube import astro
    #t0 = astro.Time()
    #day0 = int(d['mjd_day'][0])
    #sec0 = int(d['mjd_sec'][0])
    #ns0  = float(d['mjd_ns'][0])
    #t0.SetTime(day0, sec0, ns0)
    #t0.SetTime(d['mjd'][0])

    #t1 = astro.Time()
    #day1 = int(d['mjd_day'][-1])
    #sec1 = int(d['mjd_sec'][-1])
    #ns1  = float(d['mjd_ns'][-1])
    #t1.SetTime(day1, sec1, ns1)
    #t1.SetTime(d['mjd'][-1])
    #duration = t1.GetMJD() - t0.GetMJD()

    duration = d['mjd'][-1] - d['mjd'][0]
    duration *= 24 * 3600.
    return duration
