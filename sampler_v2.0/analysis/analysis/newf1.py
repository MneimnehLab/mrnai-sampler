#!/usr/bin/python

from __future__ import division
import numpy as np
import itertools
import copy
import time
import datetime
import math
import sys
import random
from collections import defaultdict
import scipy.linalg as salg
import scipy.stats as stats
import argparse
from os import listdir
from os.path import isfile, join
import sorter3
import os
from subprocess import Popen, PIPE
import single_intervals_sums_ORIG
from struct_comp import structTypeFuncs



sequences = []
# legend = {"A":"first, middle", "B":"first split, middle", "C":"first,  middle split",\
          # "D":"first, middle, last", "E":"first split, middle, last", "F":"first,  middle split, last"}






def f1scoreFromStructs(knownStruct, predictedStruct):
    # print 'computing f1 of', knownStruct, predictedStruct
    # print 'base pairs of known ='
    correct = getBasepairs(knownStruct)
    # print '=>', sorted(correct)

    #predicted structure
    # print 'base pairs of predicted ='
    predicted = getBasepairs(predictedStruct)   
    # print '=>', sorted(predicted)

    intSize = len(correct.intersection(predicted))

    precision = intSize / len(correct)
    recall = intSize / len(predicted)

    # print ''
    # print ''
    return f1score(precision, recall)


def f1score(p, r):
    return 2 * p * r / (p + r)


def alignBases(win, seq1, seq2, needToReverse=True):
    (_, i, j, w1, w2) = win
    oSeq1, oSeq2 = seq1, seq2
    # if needToReverse:
    #     seq2 = seq2[::-1]

    seq1 = "." + seq1[i-w1-1:i]
    seq2 = "." + seq2[j-w2-1:j]

    # print "finding LCS of ", seq1, seq2

    allowedPairs = [('A','U'), ('U','A'), ('C','G'), ('G','C'), ('G','U'), ('U','G')]
    m, n = len(seq1), len(seq2)
    # print 'm,n =', m,n
    H = defaultdict(int)
    backtrack = {}
    for i in range(1,m):
        for j in range(1,n):
            # print i,j, seq1[i],seq2[j]
            # if seq1[i] == seq2[j]:
            if (seq1[i], seq2[j]) in allowedPairs:
                H[(i,j)] = H[(i-1,j-1)] + 1
                backtrack[(i,j)] = (i-1,j-1)
            elif H[(i,j-1)] > H[(i-1,j)]:
                backtrack[(i,j)] = (i,j-1)
                H[(i,j)] = H[(i,j-1)]
            else:
                backtrack[(i,j)] = (i-1,j)
                H[(i,j)] = H[(i-1,j)]
            

    # print H[(m-1,n-1)]
    i = m-1
    j = n-1
    bonds = []
    while i > 0 and j > 0:
        t = backtrack[(i,j)]
        # print 't =', t
        if i==t[0]+1 and j==t[1]+1:
            # print i,j
            # re-add original indices
            bonds.append( (win[1]-win[3]-1+i, win[2]-win[4]-1+j )  )
        i = t[0]
        j = t[1]




    return bonds






def getBasepairs(S):
    global sequences

    NUM_EVEN = 2
    NUM_ODD = 2
    # print sequences
    # print 'sqs =', sequences
    # print 'getting base pairs of S =', S
    basepairs = []
    for (l,i,j,w1,w2) in S:    
        # print (l,i,j,w1,w2) 
        # print sequences[l]
        # print sequences[l+1]
        even = int(l / NUM_ODD)
        odd = l % NUM_ODD

        # print sequences[even], sequences[odd+even]
        # print (l,i,j,w1,w2)
        bonds = alignBases( (l,i,j,w1,w2)  ,  sequences[even], sequences[odd + NUM_EVEN])
        for (i_, j_) in bonds:
            basepairs.append( (l,i_,j_) )

        # print ""

    return set(basepairs)

def getBasepairIntersectionSize(S1, S2):
    s1 = getBasepairs(S1)
    s2 = getBasepairs(S2)

    return len(s1.intersection(s2))

'''
def getBasepairIntersectionSize(S1, S2):
    s1 = getBasepairs(S1)
    s2 = getBasepairs(S2)

    return len(s1.intersection(s2))
'''



def minArgMin(mylist):
    (minVal, argMin) = min( zip(mylist, range(len(mylist))), key=lambda x: x[0])
    return (minVal, argMin)


def getDistance(A,B):

    return single_intervals_sums_ORIG.processByStruct([A,B])
    


def getPathsOfFiles(parentFolder, upto=1000):
    files = []
    folder = os.path.join(parentFolder,"final_logs")
    # folder = os.path.join(parentFolder,"blah")
    # print "folder =", folder
    fnum = 0
    for qq, f in enumerate(os.listdir(folder) ):
        # print 'f =', f
        fpath = os.path.join(folder, f)
        # print 'fpath =', fpath
        if os.path.isdir(fpath):
            if upto != None and fnum >= upto:
                break
            fnum += 1

            # print 'fpath** =', fpath
            filePath = ''
            for fnum2, f2 in enumerate(os.listdir(fpath) ):
                # print 'f2 =', f2
                # begin = fileTypes[fileTypeIndx]
                begin = 'metric_cut_reprs'
                if f2[:len(begin)] == begin:
                    # print os.path.join(fpath,f2)
                    # only add to list of folders if it contains the file we are looking for
                    # need this for partially copied folders
                    files.append(os.path.join(fpath,f2))
                    break

    return files


def main():
    parser = argparse.ArgumentParser()
    # parser.add_argument('-r', '--known-data', dest='known', default='known_classes_cops')
    # parser.add_argument('-r', '--known-data', dest='known', default='_tnb_sorted_cops_wide')
    # parser.add_argument('-i', '--input-folder', dest='inputFolder', default='convg_cops_e_r20')
    parser.add_argument('-i', '--input-folder', dest='inputFolder', default='.')
    parser.add_argument('-k', '--max-to-look', dest='maxToLook', default=1)
    #parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False)
    args = parser.parse_args()

    # rnaType = 'yeast'
    rnaType = 'ribo4way_pm'
    # rnaType = 'ribo4way_exact'

    if rnaType == 'cops':
        legend = {"A":"first, middle", "B":"first split, middle", "C":"first,  middle split", "D":"first, middle, last", "E":"first split, middle, last", "F":"first,  middle split, last"}
        rfrncForClass = {k:[(0, 12, 12, 11, 11), (0, 35, 35, 8, 8), (0, 47, 47, 5, 5)] for k in ['A','B','C','D','E','F']}
        sequences = ['CGGUUUAAGUGGGCCCCGGUAAUCUUUUCGUACUCGCCAAAGUUGAAGAAGAUUAUCGGGGUUUUUGCUU', 'GCCAAAUUCACCCGGGGCCAUUAGAAAAGCAUGAGCGGUUUCAACUUCUUCUAAUAGCCCCAAAAACGAA']
    elif rnaType == 'yeast':
        ll = ['I1 detached', 'I1, I2 attached']
        jj = ['both helices', '1A', '1B', '1AB_joint']

        # kk = list(itertools.product(ll, jj))
        kk = [('I1 detached', '1A') ,('I1, I2 attached', '1A') ,('I1, I2 attached', 'both helices'),('I1 detached', 'both helices'),('I1 detached', '1AB_joint') ,('I1, I2 attached', '1AB_joint')]
        legend = {a:a for a in kk}
        # legend = {('I1, I2 attached', '1AB_joint'):('I1, I2 attached', '1AB_joint')}

    elif rnaType == 'ribo4way_pm' or rnaType == 'ribo4way_exact':
        kk = ['type1', 'type2', 'type3']
        legend = {a:a for a in kk}

        rfrncForClass = {k:[(0, 24, 20, 9, 9), (1, 10, 10, 9, 9), (2, 15, 19, 4, 4), (2, 30, 32, 5, 5), (3, 6, 6, 5, 5), (3, 14, 14, 3, 3)] for k in kk}
        
        sequences = ['UUUAUCCUGACGCUCCACCGAACG','GUCAACUCGUGGUGGCUUGC','CAGUUGAGCAUGGUCCAUUACAUGGUGCGG','AAAUAGAGAAGCGAACCAGAGAAACACACGCC']



    allTbls ={}
    for lv in legend.values():
        allTbls[lv] = list([])

    for k in range(1, int(args.maxToLook)+1):
        tblForK = getTbls(maxK = k, rnaType = rnaType, inputFolder=args.inputFolder)
        for lv in legend.values():
            allTbls[lv].append(tblForK[lv])

    print 'allTbls:'
    # print allTbls
    for lv in sorted(legend.values()):
        
        print lv,'\t','\t'.join(map(str, map(lambda y: '%.2f' %(y,), allTbls[lv] )))
    


def getTbls(maxK = 1, rnaType = 'yeast', inputFolder = ''):
    global legend
    global structType
    global sequences
    

    if rnaType == 'cops':
        legend = {"A":"first, middle", "B":"first split, middle", "C":"first,  middle split", "D":"first, middle, last", "E":"first split, middle, last", "F":"first,  middle split, last"}
        rfrncForClass = {k:[(0, 12, 12, 11, 11), (0, 35, 35, 8, 8), (0, 47, 47, 5, 5)] for k in ['A','B','C','D','E','F']}
        sequences = ['CGGUUUAAGUGGGCCCCGGUAAUCUUUUCGUACUCGCCAAAGUUGAAGAAGAUUAUCGGGGUUUUUGCUU', 'GCCAAAUUCACCCGGGGCCAUUAGAAAAGCAUGAGCGGUUUCAACUUCUUCUAAUAGCCCCAAAAACGAA']
    elif rnaType == 'yeast':
        ll = ['I1 detached', 'I1, I2 attached']
        jj = ['both helices', '1A', '1B', '1AB_joint']

        # kk = list(itertools.product(ll, jj))
        kk = [('I1 detached', '1A') ,('I1, I2 attached', '1A') ,('I1, I2 attached', 'both helices'),('I1 detached', 'both helices'),('I1 detached', '1AB_joint') ,('I1, I2 attached', '1AB_joint')]
        legend = {a:a for a in kk}
        # legend = {('I1, I2 attached', '1AB_joint'):('I1, I2 attached', '1AB_joint')}

        rfrncForClass = {}
        for (a,b) in legend.keys():
            if b == '1A':
                rfrncForClass[(a,b)] = [(0, 7, 7, 2, 2), (2, 8, 8, 3, 3), (2, 12, 13, 3, 3), (1, 16, 18, 4, 4)]
            else:
                rfrncForClass[(a,b)] = [(0, 7, 7, 2, 2), (2, 8, 8, 3, 3), (2, 12, 13, 3, 3), (1, 16, 18, 4, 4), (1, 19, 23, 2, 2)]
    
        sequences = ['NNNNUGUAUGNNNN', 'NNNNACAGAGAUGAUCAGCNNNN', 'NNNNAUGAUGUGAACUAGAUUCGNNNN', 'NNNNUACUAACACCNNNN']

    elif rnaType == 'ribo4way_pm' or rnaType == 'ribo4way_exact':
        kk = ['type1', 'type2', 'type3']
        legend = {a:a for a in kk}

        # rfrncForClass = {k:[(0, 24, 20, 9, 9), (1, 10, 10, 9, 9), (2, 15, 19, 4, 4), (2, 30, 32, 5, 5), (3, 6, 6, 5, 5), (3, 14, 14, 3, 3)] for k in kk}
        rfrncForClass = {k:[(0, 20, 24, 9, 9), (1, 10, 10, 9, 9), (2, 6, 6, 5, 5), (2, 14, 14, 3, 3), (3, 19, 15, 4, 4), (3, 32, 30, 5, 5)] for k in kk}
        
        # sequences = ['UUUAUCCUGACGCUCCACCGAACG','CAGUUGAGCAUGGUCCAUUACAUGGUGCGG', 'GUCAACUCGUGGUGGCUUGC','AAAUAGAGAAGCGAACCAGAGAAACACACGCC']
        sequences =   ["GUCAACUCGUGGUGGCUUGC","AAAUAGAGAAGCGAACCAGAGAAACACACGCC","UUUAUCCUGACGCUCCACCGAACG","CAGUUGAGCAUGGUCCAUUACAUGGUGCGG"]


    structType = structTypeFuncs[rnaType]


    
    


    clusterCountPerFile = []
    
    totalFiles = 0
    allFiles = getPathsOfFiles(inputFolder, upto=300)

    fsi = defaultdict(int)

    for numm,f in enumerate(allFiles):
        num = numm+1
        totalFiles += 1

        # read the structures in this file and 
        
        thisFileTopKStructs = []

        #print totalFiles
        #print f
        
        
        c = 0
        
        with open(f) as structsFile:
            # print f
            for line in structsFile:

                # try:
                    line = line.strip()
                    
                    data = line.split("\t")
                    # print data
                    struct = eval(data[1])
                    # print struct
                    struct = tuple(sorter3.sortConfig(struct))
                    # print struct
                    # print structType(struct)
                    # print ''
                    
                    # if c < MAX_CLUSTERS_TO_LOOK:                
                    if c < maxK:
                        thisFileTopKStructs.append( struct )

                    else:
                        break

                    c += 1
                # except:
                    # print "Unexpected error:", sys.exc_info()[0]
                    # pass        
        
        # print "count =", c
        # print ""

        # print 'thisFileTopKStructs='
        # print thisFileTopKStructs
        # clusterCountPerFile.append(c)

        for sClass in legend:
            # print 'matches for', sClass
            matches = filter(lambda x: structType(x) == sClass, thisFileTopKStructs)
            # print matches

            if len(matches) > 0:
                firstRep = matches[0]

                # print 'match for', sClass, ':', firstRep
                # if sClass == 'C':
                # print firstRep

                # fsi[(sClass, num)] = f1scoreFromStructs(firstRep, rfrncForClass[sClass])
                fsi[(sClass, num)] = f1scoreFromStructs(rfrncForClass[sClass], firstRep)

    # for k,v in fsi.iteritems():
    #     print k,v


    ret = {}

    # for l,lv in legend.iteritems():
    for l in sorted(legend.keys()):
        lv = legend[l]
        lcount = 0
        ss = 0
        
        for i in range(1, len(allFiles) + 1):
            if (l, i) in fsi:
                # print (l, i), ' -> ', fsi[(l, i)]
                lcount += 1
                ss += fsi[(l, i)]

        print '{0: <27}'.format(lv),
        if lcount > 0:
            # print '%s'.ljust %(lv,), ':',
            
            print ss/lcount
            ret[lv] = ss/lcount
        else:
            ret[lv] = 0
            print 0

    return ret





if __name__ == "__main__":

    

    main()

    # print alignBases((1, 19, 23, 7, 9), 'NNNNACAGAGAUGAUCAGCNNNN', 'NNNNAUGAUGUGAACUAGAUUCGNNNN')
