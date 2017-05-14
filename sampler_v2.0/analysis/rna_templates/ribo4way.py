sequences = 'GCAAGCCACCUCGCAGUCCUAUUU&GUCAACUCGUGGUGGCUUGC&GGCGUGGUACAUUACCUGGUACGAGUUGAC&AAAUAGAGAAGCGAACCAGAGAAACACACGCC'.split('&')

def struct_type(struct):
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

def legend():

    # kk = ["type1", "type2"]
    kk = ["type"+str(t) for t in range(1,3+1)]
    # kk = [('I1 detached', '1A') ,('I1, I2 attached', '1A') ,('I1, I2 attached', 'both helices'),('I1 detached', 'both helices'),('I1 detached', '1AB_joint') ,('I1, I2 attached', '1AB_joint')]
    legend_ = {a:a for a in kk}
        
    return legend_

def rfrncForClass():
    kk = ['type1', 'type2', 'type3']
    legend = {a:a for a in kk}

    # rfrncForClass_ = {k:[(0, 20, 24, 9, 9), (1, 10, 10, 9, 9), (2, 6, 6, 5, 5), (2, 14, 14, 3, 3), (3, 19, 15, 4, 4), (3, 32, 30, 5, 5)] for k in kk}
    rfrncForClass_ = {k:[(1, 2, 10, 10, 9, 9), (3, 0, 6, 6, 5, 5), (3, 0, 14, 14, 3, 3), (1, 0, 20, 24, 9, 9), (3, 2, 19, 15, 4, 4), (3, 2, 32, 30, 5, 5)] for k in kk}

    return rfrncForClass_

