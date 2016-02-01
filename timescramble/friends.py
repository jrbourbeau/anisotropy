from glob import glob
import os, ROOT

def makeFriends(config, masterTree='master_tree', ar='a'):

    # Keys to associate as friends
    #trees = ['ShowerLLH', 'maxLLH']
    #comps = ['proton', 'iron']
    trees = ['ShowerLLH', 'maxLLH']
    comps = ['proton', 'iron']
	#friendNames = ['_'.join([tree,comp]) for tree in trees for comp in comps]
    friendNames = ['pLLH','fLLH']
    if config in ['IT59','IT73','IT81','IT81-II']:
        dstFiles=glob('/data/user/jbourbeau/anisotropy/stripped-data/'+config+'-data/2011*.root')
    if config == 'IT81-III':
        dstFiles=glob('/data/ana/CosmicRay/Anisotropy/IceTop/IT81/2013/*.root')
    dstFiles.sort()
    llhPrefix = '/data/user/jbourbeau/anisotropy/ShowerLLH/%s-data' % (config)

    for dstFile in dstFiles:

        print 'Working on', dstFile
        if ar == 'a':
            llhFile = '/'.join([llhPrefix, os.path.basename(dstFile)])
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

    # configList = ['IT59','IT73','IT81','IT81-II','IT81-III']
    configList = ['IT73']
    addremove = raw_input('Add or remove friends [a|r]: ')
    for config in configList:
        makeFriends(config, ar=addremove)
