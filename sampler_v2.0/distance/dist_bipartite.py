#!/usr/bin/python

from __future__ import division
import numpy as np
import sys
import itertools
import shelve

MAX_DIST = 1

def winToIntervals(win):
    (l1,l2,i,j,w1,w2) = win
    return [(l1,i,w1), (l2,j,w2)]

def configToIntervals(config):
    x = []
    for c in config:
        x.extend(winToIntervals(c))
    
    return x


def process(stream):  
    
    windows = []
    configs = []

    numLines = 0
    gg = 0
    for line in stream:
        
        try:
            sp = line[:-1].split("\t")
            
            # handle both file formats
            if len(sp) >= 3:
                num, wlist, energy = sp[:3]
            else:
                wlist = line.strip()
                
            wlist = eval(wlist)
                        
            windows.extend(wlist)
            configs.append(wlist)
                
            numLines += 1
            
        except:
            pass
    
    
    matrix = np.zeros( (len(configs), len(configs)) )

    for i,j in itertools.combinations(range(len(configs)), 2):
        
        if len(configs[i]) != len(configs[j]):
            matrix[i][j] = 1 
            matrix[j][i] = 1
            
            continue
        
        
        config1 = configToIntervals(configs[i])
        config2 = configToIntervals(configs[j])
        
        intersections = 0
        unions = 0
        levelMismatch = False

        for n,(int1,int2) in enumerate(zip(config1,config2)):
            
            (lA, iA, wA) = int1
            (lB, iB, wB) = int2
    
            if lA != lB:
                levelMismatch = True
                break
                
            setA = set(range(iA-wA,iA+1))
            setB = set(range(iB-wB,iB+1))
            intersections += len(setA.intersection(setB))
            unions += len(setA.union(setB))
        
                
        
        if not levelMismatch:
            dist = 1 - (intersections / unions)
        else:
            dist = 1
        
        matrix[i][j] = dist
        matrix[j][i] = dist

    #print "XXX 4"
    
        
    #shelf = shelve.open("hier")
    #shelf["matrix"] = matrix
    #shelf.close()
    
    #print "XXX 5"
    
    NN = len(configs)
    
    for i in range(NN):
        for j in range(NN):
            print matrix[i][j],
        print ""

    
        
  
if __name__ == "__main__":
    if len(sys.argv) > 1:
        process(open(sys.argv[1]))
    else:
        process(sys.stdin)
