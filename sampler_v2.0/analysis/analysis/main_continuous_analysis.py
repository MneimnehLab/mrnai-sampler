#!/usr/bin/python

import sys
import os
from collections import defaultdict
from itertools import combinations
import copy
import compare_tbls
import numpy as np
import matplotlib.pyplot as plt
from struct_comp import structTypeFuncs
import itertools

# legend = {"A":"first, middle", "B":"first split, middle", "C":"first,  middle split", "D":"first, middle, last", "E":"first split, middle, last", "F":"first,  middle split, last"}

# appearTable = defaultdict(int)
# rankTable = defaultdict(list)


fileTypes = ["xluster", "metric_silh_reprs", "representatives"]

structType = None

try:
    parentFolder = sys.argv[1]
except:
    parentFolder = '.'

try:
    fileTypeIndx = int(sys.argv[2])
except:
    fileTypeIndx = 2

verbose = True






def doConvg(percentTablesByFolder, k):
    global structType
    global verbose

    winwMaxes = []
    n = len(percentTablesByFolder)
    # k = 25  # window size
    for i in range(1,n+1):
        end = i
        start = max(1, i-k+1)

        window = range(start, end+1)
        # print window

        mxInWin = 0
        for s,t in combinations(window, 2):
            # print s,t
            d,m = compare_tbls.findDiff(percentTablesByFolder[s-1], percentTablesByFolder[t-1])

            # x = mkTbl(1, t)
            # y = mkTbl(1, s)
            # # print x
            # # print y

            # d,m = compare_tbls.findDiff(x,y)
            # # print 'diff =',
            # # print d
            # print 'max diff=',
            # print m
            mxInWin = max(mxInWin, m)

        # print 'iter: %d,  MAX in WIN %s = %f' % (i+1, str(window), mxInWin)
        if verbose:
            print 'iter: %d,  MAX in WIN ending at %s = %f' % (i, str(window[-1]), mxInWin)
        if len(window) > 1:
            winwMaxes.append(mxInWin)

    # print 
    if verbose:
        print '\n\n'

    if verbose:
        for x in winwMaxes:
            print x

    x = range(2,n+1)
    y = winwMaxes
    plt.plot(x,y)




def getFolders(upto=None):
    folders = []
    folder = os.path.join(parentFolder,"final_logs")
    # print "folder =", folder
    for fnum, f in enumerate(os.listdir(folder) ):
        if upto != None and fnum >= upto:
            break
        fpath = os.path.join(folder, f)
        if os.path.isdir(fpath):
            filePath = ''
            for fnum2, f2 in enumerate(os.listdir(fpath) ):
                begin = fileTypes[fileTypeIndx]
                if f2[:len(begin)] == begin:

                    # only add to list of folders if it contains the file we are looking for
                    # need this for partially copied folders

                    folders.append(fpath)
                    break
            else:
                if verbose:
                    print 'X', fpath
            
    
    # for p in folders:
    #     print p

    # print 'folders len =', len(folders)


    # sys.exit(1)

    return folders[:]



# def getFolders(upto=None):
#     folders = []
#     folder = os.path.join("final_logs")
#     # print "folder =", folder
#     for fnum, f in enumerate(os.listdir(folder) ):
#         if upto != None and fnum >= upto:
#             break
#         fpath = os.path.join(folder, f)
#         if os.path.isdir(fpath):
#             filePath = ''
#             for fnum2, f2 in enumerate(os.listdir(fpath) ):
#                 begin = fileTypes[fileTypeIndx]
#                 if f2[:len(begin)] == begin:

#                     # only add to list of folders if it contains the file we are looking for
#                     # need this for partially copied folders
#                     folders.append(fpath)
#                     break
            
    
    # for p in folders:
    #     print p
    # print "hello"
    # sys.exit(1)
    

    return folders[:]


def appearsCountsForFile(folder):
    global structType
    global legend
    fpath = folder

    typeCountByNumClus = []
    typeRankByNumClus = []
    
    for numClus in range(1,10+1):
        seenTypes = []
        stRanks = {}
        typeCount = defaultdict(int)
        filePath = ''
        for fnum2, f2 in enumerate(os.listdir(fpath) ):
            fpath2 = os.path.join(fpath, f2)

            begin = fileTypes[fileTypeIndx]
            if f2[:len(begin)] == begin:
                filePath = fpath2
                break


        f = open(filePath)

        c = 0
        for line in f:

            if line[:3] == 'Opt':
                continue

            c += 1
            if c > numClus:
                break

            line = line.strip()
            data = line.split("\t")
            struct = eval(data[1])

            st = structType(struct)

            if st == None:
                continue

            if st in seenTypes:
                continue

            seenTypes.append(st)
            # print rank, data, st
            typeCount[st] += 1
            stRanks[st] = c
            
            

        f.close()

        # print typeCount
        typeCountByNumClus.append(typeCount)
        typeRankByNumClus.append(stRanks)

    # we should return transpose, so that the vals are in order: table[type][clusNums]

    table = defaultdict(list)
    rankTable = defaultdict(list)

    for numClus in range(0,10):
        x = typeCountByNumClus[numClus]
        # for k,v in x.iteritems():
        for t in legend:
            if t in x:
                v = x[t]
                table[t].append(v)
            else:
                table[t].append(0)

    for numClus in range(0,10):
        x = typeRankByNumClus[numClus]
        # for k,v in x.iteritems():
        for t in legend:
            if t in x:
                v = x[t]
                rankTable[t].append(v)
            else:
                rankTable[t].append(0)

    # print 'A'
    # print typeCountByNumClus
    # print 'B'
    # print table
    # print 'C'


    # return typeCountByNumClus
    # print 'rank table = ', rankTable
    return table, rankTable
    

def getCummulativeCounts():
    global structType
    global legend
    global verbose
    # cummulativeCounts = defaultdict(list)
    cummulativeCounts = {}
    cummulativeRanks = {}
    cummulativeRanksCount = {}
    
    percentTablesByFolder = []
    rankTablesByFolder = []

    for l in legend:
        cummulativeCounts[l] = [0]*10
        cummulativeRanks[l] = [0]*10
        cummulativeRanksCount[l] = [0]*10

    for c,fpath in enumerate(getFolders(upto=1000)):
        N = c + 1
        if N > 1:
            oldPercentages = copy.deepcopy(cummulativeCounts)

        d, ranks = appearsCountsForFile(fpath)
        for typ,countsByClus in d.iteritems():
            for i,v in enumerate(countsByClus):
                cummulativeCounts[typ][i] += v


        # print 'ranks =', ranks
        for typ,ranks in ranks.iteritems():
            for i,v in enumerate(ranks):
                # print 'xx =', typ, i, v
                cummulativeRanks[typ][i] += v
                if v > 0:
                    cummulativeRanksCount[typ][i] += 1



        # print ranks
        
        # print cummulativeCounts
        # print d
        currPercentages = {}
        currRanks = {}
        if verbose:
            print 'in dir', fpath
            print 'N =', N
        
        for k in legend.keys():
            v = cummulativeCounts[k]
            if verbose:
                print k,map(lambda y: '%.2f' %(y,), map(lambda x: float(x)/N*100, v) )
            ll = map(lambda x: float(x)/N*100, v)
            currPercentages[k] = ll
        percentTablesByFolder.append(currPercentages)

        
        for k in legend.keys():
            v = cummulativeRanks[k]
            cc = cummulativeRanksCount[k]
            # print 'v =', v
            # print 'cc =', cummulativeRanksCount[k]
            # print k,map(lambda y: '%.2f' %(y,), map(lambda x: float(x)/N*100, v) )
            ll = []
            for j in range(len(v)):
                if cc[j] > 0:
                    ll.append(float(v[j])/cc[j])
                else:
                    ll.append(0)
            currRanks[k] = ll
        rankTablesByFolder.append(currRanks)
        
        

        # if N > 1:
        #     d,m = compare_tbls.findDiff(percentTablesByFolder[-1],percentTablesByFolder[-2])
        #     print m


        # return

    # for tbl in percentTablesByFolder:
    #     for k in ['A','B','C','D','E','F']:
    #         v = tbl[k]
    #         print k,map(lambda y: '%.2f' %(y,), map(lambda x: float(x)/N*100, v) )

    #     print ''

    return percentTablesByFolder, rankTablesByFolder


def printPercentTbl(tbl):
    global legend
    # for k in ['A','B','C','D','E','F']:
    for k in sorted(legend.keys()):
        v = tbl[k]
        print legend[k],'\t','\t'.join(map(str, map(lambda y: '%.1f' %(y,), v )))
        
        
def main():
    # getTbls('yeast')
    getTbls('ribo4way_pm')
    # getTbls('ribo4way_exact')

def getTbls(rnaType = 'yeast'):
    if len(sys.argv) < 2:
        sys.stderr.write("Pass parent folder as an argument!!\n")
        return

    global legend
    global structType
    

    if rnaType == 'cops':
        legend = {"A":"first, middle", "B":"first split, middle", "C":"first,  middle split", "D":"first, middle, last", "E":"first split, middle, last", "F":"first,  middle split, last"}
    
    elif rnaType == 'yeast':
        ll = ['I1 detached', 'I1, I2 attached']
        jj = ['both helices', '1A', '1B', '1AB_joint']
        kk = list(itertools.product(ll, jj))
        # kk = [('I1 detached', '1A') ,('I1, I2 attached', '1A') ,('I1, I2 attached', 'both helices'),('I1 detached', 'both helices'),('I1 detached', '1AB_joint') ,('I1, I2 attached', '1AB_joint')]
        legend = {a:a for a in kk}

    elif rnaType == 'ribo4way':

        # kk = ["type1", "type2"]
        kk = ["type"+str(t) for t in range(1,2+1)]
        # kk = [('I1 detached', '1A') ,('I1, I2 attached', '1A') ,('I1, I2 attached', 'both helices'),('I1 detached', 'both helices'),('I1 detached', '1AB_joint') ,('I1, I2 attached', '1AB_joint')]
        legend = {a:a for a in kk}


    elif rnaType == 'ribo4way_exact' or rnaType == 'ribo4way_pm':
        
        # kk = ["type1", "type2"]
        kk = ["type"+str(t) for t in range(1,3+1)]
        # kk = [('I1 detached', '1A') ,('I1, I2 attached', '1A') ,('I1, I2 attached', 'both helices'),('I1 detached', 'both helices'),('I1 detached', '1AB_joint') ,('I1, I2 attached', '1AB_joint')]
        legend = {a:a for a in kk}
        

    structType = structTypeFuncs[rnaType]
    winSize = 25
    percentTablesByFolder, rankTablesByFolder = getCummulativeCounts()
    # return
    doConvg(percentTablesByFolder, winSize)

    print 'percent:'
    printPercentTbl(percentTablesByFolder[-1])
    print 'avg rank:'
    printPercentTbl(rankTablesByFolder[-1])

    # x = percentTablesByFolder[-1]
    # y = rankTablesByFolder[-1]
    # z = {}
    # print ''
    # print ''
    # for i in x.keys():
    #     # print x[i]
    #     # print y[i]
    #     # print ''
    #     zi = []
    #     for j in range(len(x[i])):
    #         zi.append(x[i][j])
    #         zi.append(y[i][j])
    #     z[i] = zi

    # printPercentTbl(z)


    

if __name__ == "__main__":
    main()
    