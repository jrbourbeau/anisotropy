#!/usr/bin/env python

from icecube import astro
import csv, sys

""" Simpler version calculates time difference between 1st and last event"""
def duration(d):

    t0 = astro.Time()
    day0 = int(d['mjd_day'][0])
    sec0 = int(d['mjd_sec'][0])
    ns0  = float(d['mjd_ns'][0])
    t0.SetTime(day0, sec0, ns0)

    t1 = astro.Time()
    day1 = int(d['mjd_day'][-1])
    sec1 = int(d['mjd_sec'][-1])
    ns1  = float(d['mjd_ns'][-1])
    t1.SetTime(day1, sec1, ns1)

    duration = t1.GetMJD() - t0.GetMJD()
    duration *= 24 * 3600.
    return duration


def duration_old(config, yyyymm):

    dateList = []
    t = 0.0

    if config == 'IT59':
        file = '/net/user/fmcnally/ShowerLLH/resources/IC59_GoodRunList.csv'
        f = csv.reader(open(file, 'r'))
        start, end = int(start), int(end)

        for row in f:
            date = row[0]
            if date not in ['']:
                mm = date[:2]
                dd = date[3:5]
                yy = '20'+date[6:8]
                yymmdd = int(yy+mm+dd)
                if yymmdd >= start and yymmdd <= end:
                    if row[2] != '':
                        t += float(row[2])
                    dateList.append(str(yymmdd))
        print t*60, "(s)"


    elif config == 'IT73':

        file = '/net/user/fmcnally/ShowerLLH/resources/IC79_GoodRunList.txt'
        f = open(file, 'r')
        lines = f.readlines()
        f.close()

        month = '%s-%s' % (yyyymm[:4], yyyymm[4:6])

        endList = [month in line for line in lines]
        endLine = endList.index(True)
        lines.reverse()
        startList = [month in line for line in lines]
        startLine = len(lines) - startList.index(True)
        lines.reverse()

        for line in lines[endLine:startLine]:
            t0 = line[22:26]
            IT = line[29:30]
            GOOD = line[40:44]

            if IT == 'T' and GOOD == 'GOOD':
                t += float(t0)

        return t*3600



if __name__ == "__main__":

    if len(sys.argv) != 3:
        print "Incorrect args ([config] [start yyyymm])"

    config = sys.argv[1]
    yyyymm = sys.argv[2]

    print duration_old(config, month)
