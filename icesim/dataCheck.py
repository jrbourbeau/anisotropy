#!/usr/bin/env python 

from icecube import icetray, dataio
from I3Tray import *

import sys, time, glob, os
import numpy as np


def fileCheck(fileName):

    t0 = time.time()
    tray = I3Tray()
    tray.AddModule('I3Reader', FileName=fileName)
    tray.Execute()
    tray.Finish()

    print "Time taken: ", time.time() - t0


def checker(config, out, fileList):

    badList = []
    for file in fileList:
        try: fileCheck(file)
        except RuntimeError:
            print 'Bad run found:', os.path.basename(file)
            badList += [file+'\n']

    with open(out, 'w') as f:
        f.writelines(badList)



if __name__ == "__main__":

    config = sys.argv[1]
    out = sys.argv[2]
    fileList = sys.argv[3:]

    checker(config, out, fileList)





