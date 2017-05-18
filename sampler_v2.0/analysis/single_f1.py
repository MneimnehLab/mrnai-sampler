#!/usr/bin/python

from __future__ import division
import itertools
import copy
import time
import datetime
import math
import sys
import random
from collections import defaultdict
import argparse
from os import listdir
from os.path import isfile, join
import os
# from subprocess import Popen, PIPE
# import single_intervals_sums_ORIG
# from struct_comp import structTypeFuncs
import imp



sequences = []
rfrncForClass = {}
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
    f1 = f1score(precision, recall)
    # print f1
    return f1


def f1score(p, r):
    return 2 * p * r / (p + r)


def alignBases(win, seq1, seq2, needToReverse=True):
    (l1,l2, i, j, w1, w2) = win
    oSeq1, oSeq2 = seq1, seq2
    # if needToReverse:
    #     seq2 = seq2[::-1]

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
            bonds.append( (win[2]-win[4]-1+i, win[3]-win[5]-1+j )  )
        i = t[0]
        j = t[1]



    # print bonds
    return bonds






def getBasepairs(S):
    global sequences

    NUM_EVEN = 2
    NUM_ODD = 2
    # print sequences
    # print 'sqs =', sequences
    # print 'getting base pairs of S =', S
    basepairs = []
    for (l1,l2,i,j,w1,w2) in S:    
        # print (l,i,j,w1,w2) 
        # print sequences[l]
        # print sequences[l+1]
        even = l1
        odd = l2

        # print sequences[even], sequences[odd+even]
        # print (l,i,j,w1,w2)
        bonds = alignBases( (l1,l2,i,j,w1,w2)  ,  sequences[even], sequences[odd])
        for (i_, j_) in bonds:
            basepairs.append( (l1,l2,i_,j_) )

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




def main():
    parser = argparse.ArgumentParser()
    # parser.add_argument('-r', '--known-data', dest='known', default='known_classes_cops')
    # parser.add_argument('-r', '--known-data', dest='known', default='_tnb_sorted_cops_wide')
    # parser.add_argument('-i', '--input-folder', dest='inputFolder', default='convg_cops_e_r20')
    parser.add_argument('-i', '--input-folder', dest='inputFolder', default='.')
    parser.add_argument('-k', '--max-to-look', dest='maxToLook', default=1)
    parser.add_argument('-n', '--name', dest='name', required=True)
    #parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False)
    args = parser.parse_args()

    global legend
    global structType
    global sequences
    global rfrncForClass

    rna_template_name = 'rna_templates/' + str(args.name) + '.py'
    # rna_template_name = 'rna_templates/yeast_trunc.py'
    rna_template = foo = imp.load_source('rna.template.mod', rna_template_name)

    
    structType = rna_template.struct_type
    legend = rna_template.legend()
    sequences = rna_template.sequences
    rfrncForClass = rna_template.rfrncForClass()

    
    
    f1 = f1scoreFromStructs(rfrncForClass[sClass], firstRep)




if __name__ == "__main__":

    

    main()

    # print alignBases((1, 19, 23, 7, 9), 'NNNNACAGAGAUGAUCAGCNNNN', 'NNNNAUGAUGUGAACUAGAUUCGNNNN')
