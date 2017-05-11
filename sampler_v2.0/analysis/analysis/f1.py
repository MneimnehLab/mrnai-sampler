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


def f1score(p, r):
    return 2 * p * r / (p + r)


sequences = []

def alignBases(win, seq1, seq2, needToReverse=True):
    (_, i, j, w1, w2) = win
    oSeq1, oSeq2 = seq1, seq2
    if needToReverse:
        seq2 = seq2[::-1]

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

    method = 2  # 1: without alignment.  2: with alignment

    # print 'sqs =', sequences
    basepairs = []
    for (l,i,j,w1,w2) in S:

        

        if method == 1:
            if w1 != w2:
                sys.stderr.write("unequal sized windows!! >:(  \n")
                continue

            for k in range(0,w1+1):
                i_ = i - k
                j_ = j - k
                basepairs.append( (l,i_,j_) )

        else:
            bonds = alignBases( (l,i,j,w1,w2)  ,  sequences[l], sequences[l+1])
            for (i_, j_) in bonds:
                basepairs.append( (l,i_,j_) )

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
    '''
    A = ((0, 13, 13, 12, 12), (0, 49, 49, 23, 23), (0, 70, 70, 7, 7))
    B = ((0, 4, 4, 2, 2), (0, 13, 13, 7, 7), (0, 49, 49, 23, 23))
    print "dist = ", getDistance(A,B)

    
    # sum_inters = 3 + 0 + 0
    # sum_unions = 13 + 24 + 8 + 8 + 24


    return
    '''


    parser = argparse.ArgumentParser()
    parser.add_argument('-k', '--known-data', dest='known', default='known_classes_cops')
    # parser.add_argument('-k', '--known-data', dest='known', default='_tnb_sorted_cops_wide')
    parser.add_argument('-i', '--input-folder', dest='inputFolder', default='convg_cops_e_r20')
    parser.add_argument('-m', '--max-to-look', dest='maxToLook', default=1)
    #parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False)
    args = parser.parse_args()

    MAX_CLUSTERS_TO_LOOK = int(args.maxToLook)

    knownStructures = []    #defaultdict(int)

    # read structures from known data
    with open(args.known) as knownData:
        knownStructures = [tuple(eval(line.strip())) for line in knownData]

    exactMatchCount = dict( [(struct, 0) for struct in knownStructures] )
    closeMatchDists = {}
    for struct in knownStructures:
        closeMatchDists[struct] = []
    
    
    filesOfType = 0
    

    
    '''
    for k in knownStructures:
        print k
    
    print ""
    print ""
    '''
    
    #now read files from final_folder

    #return

    clusterCountPerFile = []
    
    totalFiles = 0



    allFiles = getPathsOfFiles(args.inputFolder, upto=200)
    # for f in allFiles:
    #     print f

    # return
    



    # for f in listdir(args.inputFolder):
    #     if isfile(join(args.inputFolder,f)) and f[-5:] == ".rout":

    for numm,f in enumerate(allFiles):
            num = numm+1
            # parts = f[:-5].split("_")
            # sampler = parts[1]
            # dist = parts[2]

            #print parts
            
            # if sampler != args.sampleAlgo or dist != args.distribution:
            #     continue

            # num = int(parts[4])
            # num = int(parts[-1])

            #if  num < 360 or num > 366:
            #    continue
            totalFiles += 1

            # read the structures in this file and 
            
            thisFileStructs = []

            #print totalFiles
            #print f
            
            
            c = 0
            
            with open(f) as structsFile:
                # print f
                for line in structsFile:

                    # try:
                        line = line.strip()
                        # print line

                        
                        


                        # data = line.split(" ")
                        # #print data
                        # data = " ".join(data[1:])
                        # #print data
                        # struct = eval(data)[0]
                        
                        data = line.split("\t")
                        # print data
                        struct = eval(data[1])
                        # print struct
                        struct = tuple(sorter3.sortConfig(struct))
                        # print struct
                        # print ''
                        
                        if c < MAX_CLUSTERS_TO_LOOK:                
                            thisFileStructs.append( struct )
                        
                        #print "struct =", struct
                        c += 1
                    # except:
                        # print "Unexpected error:", sys.exc_info()[0]
                        # pass

                    #if c == MAX_CLUSTERS_TO_LOOK:
                    #    break
            

            
            # print "count =", c
            # print ""

            # print thisFileStructs
            clusterCountPerFile.append(c)

            #print "checking file", f
            # for each known structure, check if it is in this sample
            # if not, get the closest one 

            for knownStruct in knownStructures:
                # print "knownStruct =", knownStruct
                if knownStruct in thisFileStructs:
                    exactMatchCount[knownStruct] += 1
                    # print f, 'has exact struct'
                else:
                    #no exact match, so find the closest struct in this sample
                    
                    #print "not in known"
                    #ok now see which structure this is closest to


                    # for struct in thisFileStructs:
                    #     print 'struct =', struct
                    #     print 'knownStruct =', knownStruct
                    #     print 'getDist =', getDistance(knownStruct,struct) 
                    distances = [getDistance(knownStruct,struct) for struct in thisFileStructs]
                    # print 'distances =', distances
                    #print "|dist| =", len(distances)
                    if len(distances) == 0:
                        continue
                    
                    
                    (minVal, argMin) = minArgMin(distances)


                    if minVal == 1:
                        continue

                    '''
                    print "::",struct
                    for i in range(len(distances)):
                        print "X =", knownStructures[i], "dist =", distances[i]
                    '''
                    closestStruct = thisFileStructs[argMin]

                    # now check if closestStruct's closest is also knownStruct
                    
                    recipDistances = [getDistance(ks,closestStruct) for ks in knownStructures]
                    if len(recipDistances) == 0:
                        continue
                    
                    (rminVal, rargMin) = minArgMin(recipDistances)
                    closestKnownToClosest = knownStructures[rargMin]

                    '''
                    if knownStruct == ((0, 13, 13, 12, 12), (0, 49, 49, 23, 23), (0, 70, 70, 7, 7)):
                        print closestStruct
                        print "\nall structs:"
                        for struct in thisFileStructs:
                            print struct
                    '''


                    if closestKnownToClosest == knownStruct:
                        
                        '''
                        compute precision and recall
                        p = correctly predicted / correct
                        r = correctly predicted / total predicted
                        '''

                        #correct structure
                        correct = getBasepairs(knownStruct)

                        #predicted structure
                        predicted = getBasepairs(closestStruct)   

                        intSize = len(correct.intersection(predicted))

                        precision = intSize / len(correct)
                        recall = intSize / len(predicted)

                        closeMatchDists[knownStruct].append( (closestStruct, minVal, precision, recall) )


    # return

    print "totalFiles = %s" % (totalFiles,)
    

    

    # print "\t".join(["Candidate", "# hits", "avg dist", "avg f1-score", "# miss"])
    print "\t".join(["C", "% hit", "avg d", "avg f1", "% miss"])

    for n,struct in enumerate(knownStructures):

        # print "Rep", n+1
        # print "Exact:", exactMatchCount[struct]

        out = "%s\t%.3f\t" % (n+1, exactMatchCount[struct]/totalFiles)
        #\t%s\t%s\t%s
        misses = totalFiles - len(closeMatchDists[struct]) - exactMatchCount[struct]
        # print 'zxx =', len(closeMatchDists[struct])

        if len(closeMatchDists[struct]) > 0:
            avgDist      = sum([x[1] for x in closeMatchDists[struct]]) / len(closeMatchDists[struct])
            avgPrecision = sum([x[2] for x in closeMatchDists[struct]]) / len(closeMatchDists[struct])
            avgRecall    = sum([x[3] for x in closeMatchDists[struct]]) / len(closeMatchDists[struct])

            avgF1 = sum([f1score(x[2],x[3]) for x in closeMatchDists[struct]]) / len(closeMatchDists[struct])
            
            
            #misses = 100 - len(closeMatchDists[struct]) - exactMatchCount[struct]

            #print "avgDist =", avgDist
            #print "avgPrecision =", avgPrecision
            #print "avgRecall =", avgRecall

            # out += "%.3f\t%.3f" % (avgDist, f1score(avgPrecision, avgRecall) ) 
            out += "%.3f\t%.3f" % (avgDist, avgF1 ) 
            
            #print "Close: %s (%s structres) prec: %s, recall: %s, f1: %s" \
            #        % (avgDist, len(closeMatchDists[struct]), avgPrecision, avgRecall, f1score(avgPrecision, avgRecall) )
        else:
            out += "N/A\t0"
            
        out += "\t%.3f" % (misses/totalFiles, ) 

        print out
        '''
        print 
        for x in closeMatchDists[struct]:
            print x
        '''
            
    
        #print ""

    print ""
    #print "clust nums: ", clusterCountPerFile
    print "avg clust num: ", sum(clusterCountPerFile)/float(len(clusterCountPerFile))


        
    


    #print len(files)


def allLCS(seq1, seq2):
    seq1 = "." + seq1
    seq2 = "." + seq2

    # print "finding LCS of ", seq1, seq2
    m, n = len(seq1), len(seq2)
    # print 'm,n =', m,n
    H = defaultdict(int)
    # G = defaultdict(list)  # list of lists of possible bonds
    G = {}
    for i in range(-1,m):
        for j in range(-1,n):
            G[(i,j)] = set([ tuple([]) ])
    backtrack = {}
    for i in range(1,m):
        for j in range(1,n):
            # print i,j
            if seq1[i] == seq2[j] and H[(i,j-1)] == H[(i-1,j)] and H[(i,j-1)] == H[(i-1,j-1)] + 1:
                # print 'three equal'
                G[(i,j)] = set()

                for aList in G[(i-1,j-1)]:
                    G[(i,j)].add( aList + ((i,j),) )
                for aList in G[(i-1,j)]:
                    G[(i,j)].add( aList )
                for aList in G[(i,j-1)]:
                    G[(i,j)].add( aList  )
                
                H[(i,j)] = H[(i-1,j-1)] + 1

            elif seq1[i] == seq2[j]:
                # print 'A'
                G[(i,j)] = set()
                H[(i,j)] = H[(i-1,j-1)] + 1
                for aList in G[(i-1,j-1)]:
                    # print 'alist =', aList
                    G[(i,j)].add( aList  + ((i,j),)  )

            elif H[(i,j-1)] == H[(i-1,j)]:
                # print 'two equal'
                G[(i,j)] = set()
                for aList in G[(i-1,j)]:
                    G[(i,j)].add( aList )
                for aList in G[(i,j-1)]:
                    G[(i,j)].add( aList)
                H[(i,j)] = H[(i,j-1)]

            elif H[(i,j-1)] > H[(i-1,j)]:
                # print 'B'
                G[(i,j)] = set()
                H[(i,j)] = H[(i,j-1)]
                for aList in G[(i,j-1)]:
                    G[(i,j)].add( aList )
            else:
                # print 'C'
                G[(i,j)] = set()
                H[(i,j)] = H[(i-1,j)]
                for aList in G[(i-1,j)]:
                    G[(i,j)].add( aList )

            print 'chose =', H[(i,j)]
            # print 'lists =', [list(x) for x in G[(i,j)] ]
            # print ''

    print 'last ='
    for x in G[(m-1,n-1)]:
        print x
    print 



            

if __name__ == "__main__":
    '''    
    S1 = [(0, 7, 7, 2, 2), (2, 8, 8, 3, 3), (2, 12, 13, 3, 3), (1, 16, 18, 4, 4), (1, 19, 23, 1, 1)]
    S2 = [(0, 7, 7, 2, 2), (2, 9, 9, 3, 3), (2, 12, 13, 2, 2), (1, 16, 18, 4, 4), (1, 19, 23, 2, 2)]
    
    print sorted(getBasepairs(S1)), len(getBasepairs(S1))
    print sorted(getBasepairs(S2)), len(getBasepairs(S2))

    print getBasepairIntersectionSize(S1, S2)
    '''

    sequences = ['CGGUUUAAGUGGGCCCCGGUAAUCUUUUCGUACUCGCCAAAGUUGAAGAAGAUUAUCGGGGUUUUUGCUU', 'AAGCAAAAACCCCGAUAAUCUUCUUCAACUUUGGCGAGUACGAAAAGAUUACCGGGGCCCACUUAAACCG']

    main()


    # LCS('ABAB', 'BABA')
    # allLCS('AAVBBCCX','CCVBBAAX')
    


    seq1 = 'CGGUUUAAGUGGGCCCCGGUAAUCUUUUCGUACUCGCCAAAGUUGAAGAAGAUUAUCGGGGUUUUUGCUU'
    seq2 = 'AAGCAAAAACCCCGAUAAUCUUCUUCAACUUUGGCGAGUACGAAAAGAUUACCGGGGCCCACUUAAACCG'



    # alignBonds((0,37,35,12,11), 'CGGUUUAAGUGGGCCCCGGUAAUCUUUUCGUACUCGCCAAAGUUGAAGAAGAUUAUCGGGGUUUUUGCUU', 'AAGCAAAAACCCCGAUAAUCUUCUUCAACUUUGGCGAGUACGAAAAGAUUACCGGGGCCCACUUAAACCG')
    # print alignBases((0,9,5,8,4), 'YCCXXXCCC', 'GGGGG')
    # print alignBases((0,37,35,12,11), 'CGGUUUAAGUGGGCCCCGGUAAUCUUUUCGUACUCGCCAAAGUUGAAGAAGAUUAUCGGGGUUUUUGCUU', 'AAGCAAAAACCCCGAUAAUCUUCUUCAACUUUGGCGAGUACGAAAAGAUUACCGGGGCCCACUUAAACCG')


