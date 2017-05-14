sequences = 'ACCAGAGAAACACACGAAAAAAAAAGAGAGAAGUGAA&CGUGGUAUAUUACCUGGU&UCACAGUCCUCUUU'.split('&')

def struct_type(struct):
    # struct = sorted(struct)

    DELTA = 3

    if len(struct) == 3:
        type1Std = [(0, 1, 5, 5, 4, 4), (0, 1, 16, 18, 3, 3), (0, 2, 36, 14, 13, 13)]

        allTrue = True
        for pWin, cWin in zip(struct, type1Std):  # pWin = predicted, cWin = correct
            if pWin[0] == cWin[0] and pWin[1] == cWin[1] and\
                pWin[2] - pWin[4] >= cWin[2] - cWin[4] - DELTA and\
                pWin[3] - pWin[5] >= cWin[3] - cWin[5] - DELTA and\
                pWin[2] <= cWin[2] + DELTA and\
                pWin[3] <= cWin[3] + DELTA:
                # print 'yes for', pWin, cWin
                pass
            else:
                # print 'no for', pWin, cWin
                allTrue = False

        if allTrue:
            return "type1"

    elif len(struct) == 4:
        # type2Std = [(0, 24, 20, 9, 9), (1, 10, 10, 9, 9), (2, 15, 19, 4, 4), (2, 30, 32, 5, 5), (3, 6, 6, 5, 5), (3, 14, 14, 3, 3)]
        # type2Std =  [(0, 5, 5, 4, 4), (0, 16, 18, 3, 3), (1, 28, 6, 5, 5), (1, 36, 14, 3, 3)]
        type2Std = [(0, 1, 5, 5, 4, 4), (0, 1, 16, 18, 3, 3), (0, 2, 28, 6, 5, 5), (0, 2, 36, 14, 3, 3)]
        
        
        allTrue = True
        for pWin, cWin in zip(struct, type2Std):  # pWin = predicted, cWin = correct
            if pWin[0] == cWin[0] and pWin[1] == cWin[1] and\
                pWin[2] - pWin[4] >= cWin[2] - cWin[4] - DELTA and\
                pWin[3] - pWin[5] >= cWin[3] - cWin[5] - DELTA and\
                pWin[2] <= cWin[2] + DELTA and\
                pWin[3] <= cWin[3] + DELTA:
                # print 'yes for', pWin, cWin
                pass
            else:
                # print 'no for', pWin, cWin
                allTrue = False

        if allTrue:
            # if struct[4][3] == struct[4][4]:
            # if struct[1][3] != struct[1][4] or struct[2][3] != struct[2][4]:
            if struct[0][4] == struct[0][5] and struct[1][4] == struct[1][5] and struct[2][4] == struct[2][5] and struct[3][4] == struct[3][5]:
                return "type2"
            else:
                # print struct
                return "type3"

    return None

def legend():
    # kk = ["type1", "type2"]
    kk = ["type"+str(t) for t in range(1,3+1)]
    # kk = [('I1 detached', '1A') ,('I1, I2 attached', '1A') ,('I1, I2 attached', 'both helices'),('I1 detached', 'both helices'),('I1 detached', '1AB_joint') ,('I1, I2 attached', '1AB_joint')]
    legend_ = {a:a for a in kk}

    return legend_

def rfrncForClass():
    kk = ["type"+str(t) for t in range(1,3+1)]
    # [(0, 1, 5, 5, 4, 4), (0, 1, 16, 18, 3, 3), (0, 2, 28, 6, 5, 5), (0, 2, 36, 14, 3, 3)]
    return {k:[(0, 1, 5, 5, 4, 4), (0, 1, 16, 18, 3, 3), (0, 2, 28, 6, 5, 5), (0, 2, 36, 14, 3, 3)] for k in kk}



