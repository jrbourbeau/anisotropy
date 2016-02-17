#!/usr/bin/env python

import argparse, os, ROOT

import myGlobals as my
import dataFunctions as df

def makeFriends(config, masterTree='master_tree', ar='a'):

    # Setup global paths
    my.setupShowerLLH(verbose=False)

    # Keys to associate as friends
    friendNames = ['ML_energy','pLLH','hLLH','oLLH','fLLH','llh_comp']

    # Collect DST and ShowerLLH root files
    dstFiles = df.getDSTfiles(config)
    llhPrefix = '%s/%s_data/ani_files' % (my.llh_data, config)

    for dstFile in dstFiles:

        print 'Working on', dstFile
        if ar == 'a':
            llhFile = '%s/%s' % (llhPrefix, os.path.basename(dstFile))
            if not os.path.isfile(llhFile):
                print llhFile, 'not found...'
                continue

        f = ROOT.TFile(dstFile, 'update')
        t = f.Get(masterTree)
        for friend in friendNames:
            # Option to remove friends instead
            if ar == 'r':
                fr = t.GetListOfFriends().FindObject(friend)
                t.GetListOfFriends().Remove(fr)
            if ar == 'a':
                t.AddFriend(friend, llhFile)
        t.Write()
        f.Close()



if __name__ == "__main__":

    p = argparse.ArgumentParser(
            description='Links ShowerLLH files to anisotropy root files')
    p.add_argument('-c', '--config', dest='config', nargs='*',
            default=['IT73','IT81','IT81-II','IT81-III','IT81-IV'],
            help='Detector configurations to run over')
    p.add_argument('--ar', dest='addremove',
            choices=['a','r'],
            help='Add (a) or remove (r) friends')
    args = p.parse_args()

    if args.addremove == None:
        raise SystemExit('Need to specify "--ar" option. Exiting...')

    for config in args.config:
        makeFriends(config, ar=args.addremove)



