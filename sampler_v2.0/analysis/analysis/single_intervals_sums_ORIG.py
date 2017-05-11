from __future__ import division
import numpy as np
import sys
import itertools
import shelve

MAX_DIST = 1

def winToIntervals(win):
    (l,i,j,w1,w2) = win
    return [(l,i,w1), (l+1,j,w2)]

def configToIntervals(config):
    x = []
    for c in config:
        x.extend(winToIntervals(c))
    
    return x

'''
FIX ASAP!!!!
'''
def processByStruct(listOfStructs):  
    windows = []
    configs = []
    kept = []
    keptCount = 0
    numLines = 0
    gg = 0
    #for line in open(fname):
        #if numLines >= 3250: break

    for struct_ in listOfStructs:
        
        try:
            
            wlist = struct_
            unequal = len([1 for l,i,j,w1,w2 in wlist if w1 != w2])
            #unequal = 0
            
            #check for gapless adjacent windows
            gapless = False
            for k in range(len(wlist)-1):
                l,i,j,w1,w2 = wlist[k]
                lx,ix,jx,w1x,w2x = wlist[k+1]
                if i + 1 == ix - w1x and j + 1 == jx - w2x:
                    gapless = True
            
            if gapless:
                #sys.stderr.write( "%s gapless! :: %s \n" % (gg, wlist ))
                gg += 1
                
            if True or unequal == 0:  # and not gapless:
                windows.extend(wlist)
                configs.append(wlist)
                kept.append(numLines)
                keptCount += 1
            else:
                return 1
                
            if keptCount >= 5000:
                break
            '''
            print num
            print wlist
            print energy
            print ""
            '''
            numLines += 1
            
            #print struct_
            
        except:
            #sys.stderr.write("error\n")
            pass

    #print ""

    windows = list(set(windows))
    windows.sort()
    
    intervals = []
    for (l, i, j, w1, w2) in windows:
        interval1 = (l,i,w1)
        interval2 = (l+1,j,w2)
        intervals.append(interval1)
        intervals.append(interval2)
        
    intervals = list(set(intervals))
    intervals.sort()
    
    
    
    
    reverseMap = {i:n  for n,i in enumerate(intervals)}
        
    similarities = np.ones( (len(intervals),len(intervals) ) )
    
    allSims = []
    for a,b in itertools.combinations(range(len(intervals)), 2):
       
        (lA, iA, wA) = intervals[a]         #(l,  10, 4) =>   [6,7,8,9,10]
        (lB, iB, wB) = intervals[b]
        
        if lA != lB:
            jaccardSim = 0
            
        else:
            setA = set(range(iA-wA,iA+1))
            setB = set(range(iB-wB,iB+1))
            jaccardSim = len(setA.intersection(setB)) / len(setA.union(setB))
            
            
        similarities[a][b] = jaccardSim
        similarities[b][a] = jaccardSim
        
        allSims.append( (a,b,jaccardSim) )

    #print "configs =", configs
    #print (len(configs), len(configs))
    matrix = np.zeros( (len(configs), len(configs)) )

    # now look pick any two configs:  0 and 4
    for i,j in itertools.combinations(range(len(configs)), 2):
        
        if len(configs[i]) != len(configs[j]):
            matrix[i][j] = 1 
            matrix[j][i] = 1
            
            continue
        
        
        config1 = configToIntervals(configs[i])
        config2 = configToIntervals(configs[j])
            
        
        prod = 1
        ss = 0
        intersections = 0
        unions = 0
        prt = False
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
        
        
            intId1 = reverseMap[int1]
            intId2 = reverseMap[int2]
            
        
        if not levelMismatch:
            dist = 1 - (intersections / unions)
        else:
            dist = 1
        
        avg = ss/len(config1)
        
        matrix[i][j] = dist
        matrix[j][i] = dist

    
        
    #shelf = shelve.open("hier")
    #shelf["matrix"] = matrix
    #shelf.close()
    
    
    NN = len(configs)
    '''
    for i in range(NN):
        for j in range(NN):
            print matrix[i][j],
        print ""
    '''
    #print matrix


    #print ""
    return matrix[0][1]

    
'''
This is optimized to work with more data
'''
def process(fname, mapOnly = False):  
    #print "here"
    windows = []
    configs = []
    kept = []
    keptCount = 0
    numLines = 0
    gg = 0
    for line in open(fname):
        if numLines >= 7250: break
        
        try:
            sp = line[:-1].split("\t")
            
            if len(sp) >= 3:
                #wlist = eval(line.strip())
                num, wlist, energy = sp[:3]
                
            else:
                wlist = line.strip()
                
            wlist = eval(wlist)
            #num = int(num)
            #wlist = eval(wlist)
            #energy = float(energy)
            #wlist = eval(line.strip())
            
            #print wlist
            unequal = len([1 for l,i,j,w1,w2 in wlist if w1 != w2])
            #unequal = 0
            
            #check for gapless adjacent windows
            gapless = False
            for k in range(len(wlist)-1):
                l,i,j,w1,w2 = wlist[k]
                lx,ix,jx,w1x,w2x = wlist[k+1]
                if i + 1 == ix - w1x and j + 1 == jx - w2x:
                    gapless = True
            
            if gapless:
                sys.stderr.write( "%s gapless! :: %s \n" % (gg, wlist ))
                gg += 1
                
            #if unequal != 0:
            #    sys.stderr.write( "unequal! :: %s \n" % (wlist ))
                #print 
            
            #print unequal
            
            #if not gapless:
            #if unequal == 0:
            if True or (unequal == 0 and not gapless):
                windows.extend(wlist)
                configs.append(wlist)
                kept.append(numLines)
                keptCount += 1
                
            if keptCount >= 7000:
                break
            '''
            print num
            print wlist
            print energy
            print ""
            '''
            numLines += 1
            
            
            
        except:
            #sys.stderr.write("error\n")
            pass
    
    
    #for n,k in enumerate(kept):
    #    sys.stderr.write("%s > %s \n" %( n,k))

    #print "gapless =", gg
    if mapOnly:
        for n,k in enumerate(kept):
            print "%s > %s" %( n,k)

        return
    
    #print "XXX 1"
    
    windows = list(set(windows))
    windows.sort()
    '''
    for n,w in enumerate(windows):
        print n,w
    '''
    
    intervals = []
    for (l, i, j, w1, w2) in windows:
        interval1 = (l,i,w1)
        interval2 = (l+1,j,w2)
        intervals.append(interval1)
        intervals.append(interval2)
        
    intervals = list(set(intervals))
    intervals.sort()
    
    #print "XXX 2"
    
    
    reverseMap = {i:n  for n,i in enumerate(intervals)}
        
    similarities = np.ones( (len(intervals),len(intervals) ) )
    
    allSims = []
    for a,b in itertools.combinations(range(len(intervals)), 2):
       
        (lA, iA, wA) = intervals[a]         #(l,  10, 4) =>   [6,7,8,9,10]
        (lB, iB, wB) = intervals[b]
        
        if lA != lB:
            jaccardSim = 0
            
        else:
            setA = set(range(iA-wA,iA+1))
            setB = set(range(iB-wB,iB+1))
            jaccardSim = len(setA.intersection(setB)) / len(setA.union(setB))
            
            
        similarities[a][b] = jaccardSim
        similarities[b][a] = jaccardSim
        
        allSims.append( (a,b,jaccardSim) )
    
    '''
    allSims.sort( key=lambda x: x[2] )
    
    for a,b,s in allSims:
        print intervals[a], intervals[b], s
    '''
    #return
    #print "XXX 3"

    matrix = np.zeros( (len(configs), len(configs)) )

    # now look pick any two configs:  0 and 4
    for i,j in itertools.combinations(range(len(configs)), 2):
        
        if len(configs[i]) != len(configs[j]):
            matrix[i][j] = 1 
            matrix[j][i] = 1
            
            continue
        
        
        config1 = configToIntervals(configs[i])
        config2 = configToIntervals(configs[j])
            
        
        prod = 1
        ss = 0
        intersections = 0
        unions = 0
        prt = False
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
        
        
            intId1 = reverseMap[int1]
            intId2 = reverseMap[int2]
            
        
        if not levelMismatch:
            dist = 1 - (intersections / unions)
        else:
            dist = 1
        
        avg = ss/len(config1)
        
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
    if len(sys.argv) > 2:
        if "-map" in sys.argv:
            process(sys.argv[1], True)
    elif len(sys.argv) > 1:
        process(sys.argv[1])
    else:
        print "need file name!"
    