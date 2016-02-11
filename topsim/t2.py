#!/usr/bin/env python

import glob

if __name__ == "__main__":

    with open('badFiles.txt','r') as f:
        badFiles = f.readlines()
        badFiles = [file.strip() for file in badFiles]

    exeFiles = glob.glob('/home/fmcnally/npx4/npx4-execs/*.sh')
    exeFiles.sort()

    for exeFile in exeFiles:

        with open(exeFile,'r') as f:
            exeLines = f.readlines()

        ex = exeLines[-2].strip().split(' ')
        for badFile in badFiles:
            if badFile in ex:
                ex.remove(badFile)

        exeLines[-2] = ' '.join(ex) + '\n'
        with open(exeFile,'w') as f:
            f.writelines(exeLines)
