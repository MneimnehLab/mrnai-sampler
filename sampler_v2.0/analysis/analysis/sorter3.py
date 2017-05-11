import sys


'''
Cost of computing the boundary is O(|W|)?
'''
def getConfigBoundary(C, numRNAs):
    winsOnLevels = {}
    for l in range(numRNAs - 1):
        winsOnLevels[l] = []
    for w in C:
        (l,i,j,w1,w2) = w
        winsOnLevels[l].append(w)

    for l in range(numRNAs - 1):
        winsOnLevels[l].sort()

    boundary = [0] * numRNAs

    for l in range(0, numRNAs - 1):
        if winsOnLevels[l] != []:
            (l,i,j,w1,w2) = winsOnLevels[l][-1]
            boundary[l] = max(boundary[l], i)
            boundary[l+1] = max(boundary[l+1], j)

    if winsOnLevels[numRNAs - 2] != []:
        (l,i,j,w1,w2) = winsOnLevels[numRNAs - 2][-1]
        boundary[numRNAs - 1] = max(boundary[numRNAs - 1], j)

    return boundary, winsOnLevels



def getTerminalWindowFromUnsorted(S, numRNAs, boundary, winsOnLevels):
    #boundary, winsOnLevels = getConfigBoundary(S, numRNAs)

    terminal = None
    for l in range(numRNAs-2, -1, -1):

        if len(winsOnLevels[l]) > 0:
            lastWin = winsOnLevels[l][-1]
            (l,i,j,u,v) = lastWin

            if i == boundary[l] and j == boundary[l+1]:
                terminal = lastWin

                #update boundary and winsOnLevels
                #try to do this without o(|W|) cost

                break

    



    return terminal

def sortConfig(S):
    configs = [ [] ]
    numRNAs = 6
    done = False

    reverseSortedS = []
    
    #print "tosort = " ,S
    boundary, winsOnLevels = getConfigBoundary(S, numRNAs)

    C = list(S)
    while len(C) > 0:
        
        #print "C =", C

        terminalInC = getTerminalWindowFromUnsorted(C, numRNAs, boundary, winsOnLevels)
        reverseSortedS.append(terminalInC)

        C = [win for win in C if win != terminalInC]

        boundary, winsOnLevels = getConfigBoundary(C, numRNAs)


    return reverseSortedS[::-1]




def main():
    S = [ (1,30,30,3,3), (2, 12, 13, 3, 3), (1, 16, 18, 4, 4), (0, 7, 7, 2, 2), (1, 19, 23, 2, 2), (2, 8, 8, 3, 3)]
    S = [(2, 12, 13, 3, 3), (1, 16, 18, 4, 4), (0, 7, 7, 2, 2), (1, 19, 23, 2, 2), (2, 8, 8, 3, 3)]
    #S = [(2, 12, 13, 3, 3), (0, 7, 7, 2, 2)]

    #print sortConfig(S)


    for line in sys.stdin:
        data = line.strip().split("\t")
        
        if len(data) == 1:
            S = eval(data[0])
            print sortConfig(S)

        elif len(data) >= 3:
            S = eval(data[1])
            data[1] = sortConfig(S)
            print "\t".join(map(str,data))



    

if __name__ == "__main__":
    main()

