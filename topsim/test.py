#!/usr/bin/env python

from icecube import dataio
import glob

if __name__ == "__main__":

    exeFiles = glob.glob('/home/fmcnally/npx4/npx4-execs/*')
    exeFiles.sort()
    testFiles, badFiles = [],[]

    for file in exeFiles:

        with open(file, 'r') as f:
            lines = f.readlines()
        lines = lines[-2].strip().split(' ')
        lines = [l for l in lines if l[-7:]=='.i3.bz2']

        testFiles += lines
        
    for i, file in enumerate(testFiles):

        if i%100 == 0:
            print '%i files checked...' % i

        i3File = dataio.I3File(file)
        framecount = -1
        while i3File.more():
            framecount += 1
            try: frame = i3File.pop_frame()
            except Exception as e:
                badFiles += [file]
                print 'File:', file
                print 'Frame:', frame
                print 'Error:', e
                i3File.close()
                break

        

