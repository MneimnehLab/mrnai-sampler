# legend = {"A":"first, middle", "B":"first split, middle", "C":"first,  middle split",\
#           "D":"first, middle, last", "E":"first split, middle, last", "F":"first,  middle split, last"}


def structTypeCops(struct):

    MIDDLE_RIGHT_WIN_END = 50
    MIDDLE_LEFT_WIN_START = 24

    if len(struct) == 2:
        A = struct[0]
        B = struct[1]
        
        # first, middle
        if                                         A[1] <= 15 and A[2] <= 15 and\
           B[1]-B[3] >= 24 and B[2]-B[4] >= 24 and B[1] <= 50 and B[2] <= 50:
            return "A"

    elif len(struct) == 3:
        A = struct[0]
        B = struct[1]
        C = struct[2]

        # first split, middle
        if                                         B[1] <= 15 and B[2] <= 15 and    \
           C[1]-C[3] >= 24 and C[2]-C[4] >= 24 and C[1] <= 50 and C[2] <= 50:
            return "B"      
          
        
        #first,   middle split
        elif A[1] <= 15 and A[2] <= 15 and \
            B[1]-B[3] >= 25 and B[2]-B[4] >= 25 and B[1] <= MIDDLE_RIGHT_WIN_END and B[2] <= MIDDLE_RIGHT_WIN_END and\
            C[1]-C[3] >= MIDDLE_LEFT_WIN_START and C[2]-C[4] >= MIDDLE_LEFT_WIN_START and C[1] <= 50 and C[2] <= 50:
            return "C"

            
        # first, middle, last
        elif                                         A[1] <= 15 and A[2] <= 15 and\
             B[1]-B[3] >= 24 and B[2]-B[4] >= 24 and B[1] <= 50 and B[2] <= 50 and\
             C[1]-C[3] >= 60 and C[2]-C[4] >= 60 and C[1] <= 70 and C[2] <= 70:
            return "D"

    elif len(struct) == 4:
        A = struct[0]
        B = struct[1]
        C = struct[2]
        D = struct[3]

        # first split, middle, last
        if                                         B[1] <= 15 and B[2] <= 15 and \
           C[1]-C[3] >= 24 and C[2]-C[4] >= 24 and C[1] <= 50 and C[2] <= 50 and \
           D[1]-D[3] >= 60 and D[2]-D[4] >= 60 and D[1] <= 70 and D[2] <= 70:
            return "E"

        # first,  middle split, last
        if                                         A[1] <= 15 and A[2] <= 15 and \
           B[1]-B[3] >= 24 and B[2]-B[4] >= 24 and B[1] <= MIDDLE_RIGHT_WIN_END and B[2] <= MIDDLE_RIGHT_WIN_END and \
           C[1]-C[3] >= 39 and C[2]-C[4] >= 39 and C[1] <= 50 and C[2] <= 50 and \
           D[1]-D[3] >= 60 and D[2]-D[4] >= 60 and D[1] <= 70 and D[2] <= 70:
            
            return "F"



def structTypeYeast(struct):
    introns = ""
    helix = ""
    level0 = [(l,i,j,w1,w2) for (l,i,j,w1,w2) in struct if l == 0]
    level1 = [(l,i,j,w1,w2) for (l,i,j,w1,w2) in struct if l == 1]
    level2 = [(l,i,j,w1,w2) for (l,i,j,w1,w2) in struct if l == 3]

    if level0 == [] and level2 == []:
        introns = "I1, I2 detached"
    elif  level0 == [] or level0 == [(0, 10, 6, 1, 1)]:
        introns = "I1 detached"
    elif  level2 == []:
        introns = "I2 detached"
    else:
        introns = "I1, I2 attached"

    has1A = False
    has1B = False
    has1A_B = False
    for win in level1:
        if win == (1, 16, 18, 4, 4):
            has1A = True
        if win == (1, 19, 23, 2, 2):
            has1B = True 
        if win == (1, 19, 23, 7, 9):
            has1A_B = True


    if has1B and has1A:
        helix = "both helices"
    elif has1B:
        helix = "1B"
    elif has1A:
        helix = "1A"
    elif has1A_B:
        helix = "1AB_joint"
    else:
        helix = "no helix"

    return (introns, helix)
    


def structTypeRibo4Way_Exact(struct):
    # return None
    # type 1:  [(0, 24, 20, 9, 9), (1, 10, 10, 9, 9), (2, 15, 19, 4, 4), (2, 30, 32, 5, 5), (3, 14, 14, 13, 13)]
    # type 2:  [(0, 24, 20, 9, 9), (1, 10, 10, 9, 9), (2, 15, 19, 4, 4), (2, 30, 32, 5, 5), (3, 6, 6, 5, 5), (3, 14, 14, 3, 3)]
    # type 3:  [(0, 24, 20, 9, 9), (1, 10, 10, 9, 9), (2, 15, 19, 4, 4), (2, 30, 32, 5, 5), (3, 9, 8, 8, 7), (3, 14, 14, 3, 3)]


    # type1 = set([(0, 24, 20, 9, 9), (1, 10, 10, 9, 9), (2, 15, 19, 4, 4), (2, 30, 32, 5, 5), (3, 14, 14, 13, 13)])
    # type2 = set([(0, 24, 20, 9, 9), (1, 10, 10, 9, 9), (2, 15, 19, 4, 4), (2, 30, 32, 5, 5), (3, 6, 6, 5, 5), (3, 14, 14, 3, 3)])
    allTypes = [
        [(0, 24, 20, 9, 9), (1, 10, 10, 9, 9), (2, 15, 19, 4, 4), (2, 30, 32, 5, 5), (3, 14, 14, 13, 13)],
        [(0, 24, 20, 9, 9), (1, 10, 10, 9, 9), (2, 15, 19, 4, 4), (2, 30, 32, 5, 5), (3, 6, 6, 5, 5), (3, 14, 14, 3, 3)],
        [(0, 24, 20, 9, 9), (1, 10, 10, 9, 9), (2, 15, 19, 4, 4), (2, 30, 32, 5, 5), (3, 9, 8, 8, 7), (3, 14, 14, 3, 3)]
    ]

    for i,st in enumerate(allTypes):
        s = set(struct)
        if s == set(st):
            return "type"+str(i+1)

    return None



def structTypeRibo4Way_PlusMinus(struct):
    # struct = sorted(struct)

    # print "struct to compare =", struct

    if len(struct) == 5:
        # type1Std = [(0, 24, 20, 9, 9), (1, 10, 10, 9, 9), (2, 15, 19, 4, 4), (2, 30, 32, 5, 5), (3, 14, 14, 13, 13)]
        # type1Std = [(0, 20, 24, 9, 9), (1, 10, 10, 9, 9), (2, 14, 14, 13, 13), (3, 19, 15, 4, 4), (3, 32, 30, 5, 5)]
        type1Std = [(1, 2, 10, 10, 9, 9), (3, 0, 14, 14, 13, 13), (1, 0, 20, 24, 9, 9), (3, 2, 19, 15, 4, 4), (3, 2, 32, 30, 5, 5)]

        allTrue = True
        for pWin, cWin in zip(struct, type1Std):  # pWin = predicted, cWin = correct
            if pWin[0] == cWin[0] and pWin[1] == cWin[1] and\
                pWin[2] - pWin[4] >= cWin[2] - cWin[4] - 3 and\
                pWin[3] - pWin[5] >= cWin[3] - cWin[5] - 3 and\
                pWin[2] <= cWin[2] + 3 and\
                pWin[3] <= cWin[3] + 3:
                # print 'yes for', pWin, cWin
                pass
            else:
                # print 'no for', pWin, cWin
                allTrue = False

        if allTrue:
            # print 'type1'
            return "type1"

    elif len(struct) == 6:
        # type2Std = [(0, 24, 20, 9, 9), (1, 10, 10, 9, 9), (2, 15, 19, 4, 4), (2, 30, 32, 5, 5), (3, 6, 6, 5, 5), (3, 14, 14, 3, 3)]
        # type2Std = [(0, 20, 24, 9, 9), (1, 10, 10, 9, 9), (2, 6, 6, 5, 5), (2, 14, 14, 3, 3), (3, 19, 15, 4, 4), (3, 32, 30, 5, 5)]
        type2Std = [(1, 2, 10, 10, 9, 9), (3, 0, 6, 6, 5, 5), (3, 0, 14, 14, 3, 3), (1, 0, 20, 24, 9, 9), (3, 2, 19, 15, 4, 4), (3, 2, 32, 30, 5, 5)]

        allTrue = True
        for pWin, cWin in zip(struct, type2Std):  # pWin = predicted, cWin = correct
            if pWin[0] == cWin[0] and pWin[1] == cWin[1] and\
                pWin[2] - pWin[4] >= cWin[2] - cWin[4] - 3 and\
                pWin[3] - pWin[5] >= cWin[3] - cWin[5] - 3 and\
                pWin[2] <= cWin[2] + 3 and\
                pWin[3] <= cWin[3] + 3:
                # print 'yes for', pWin, cWin
                pass
            else:
                # print 'no for', pWin, cWin
                allTrue = False

        if allTrue:
            # if struct[4][3] == struct[4][4]:
            if struct[1][4] == struct[1][5]:
                # print 'type2'
                return "type2"
            else:
                # print struct
                # print 'type3'
                return "type3"

    return None

    


structTypeFuncs = {'yeast':structTypeYeast, 'cops': structTypeCops, 'ribo4way_exact': structTypeRibo4Way_Exact, 'ribo4way_pm': structTypeRibo4Way_PlusMinus}



if __name__ == '__main__':
    s = [(0, 24, 20, 9, 9), (1, 10, 10, 9, 9), (2, 15, 19, 4, 4), (2, 30, 32, 5, 5), (3, 9, 8, 8, 7), (3, 14, 14, 3, 3)]
    structTypeRibo4Way_PlusMinus(s)

