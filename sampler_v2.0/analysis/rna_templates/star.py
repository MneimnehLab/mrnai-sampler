sequences = 'UGACGGGUUUGG&ACCUUAUUUGAA&CCAAACCCGUCAAUCAAGUCUACACUGUUCAAAUAAGGU&CAGUGUACUGAUGAGUCCGUGAGGACGAAACUUGAU'.split('&')

def struct_type(struct):
    '''
    # print 's =', struct
    # if struct == [(0, 12, 12, 11, 11), (1, 39, 12, 11, 11), (2, 19, 7, 6, 6), (2, 27, 36, 6, 6)]:
    if struct == [(2, 0, 12, 12, 11, 11), (2, 3, 19, 7, 6, 6), (2, 3, 27, 36, 6, 6), (2, 1, 39, 12, 11, 11)]:
        return 'type1'
    else:
        return None
    '''
    DELTA = 3

    if len(struct) == 4:
        # type1Std = [(0, 24, 20, 9, 9), (1, 10, 10, 9, 9), (2, 15, 19, 4, 4), (2, 30, 32, 5, 5), (3, 14, 14, 13, 13)]
        # type1Std = [(0, 12, 12, 11, 11), (1, 39, 12, 11, 11), (2, 19, 7, 6, 6), (2, 27, 36, 6, 6)]
        type1Std = [(2, 0, 12, 12, 11, 11), (2, 3, 19, 7, 6, 6), (2, 3, 27, 36, 6, 6), (2, 1, 39, 12, 11, 11)]

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
            # print 'ret type1'
            return "type1"



    if len(struct) == 5:
        # type1Std = [(0, 24, 20, 9, 9), (1, 10, 10, 9, 9), (2, 15, 19, 4, 4), (2, 30, 32, 5, 5), (3, 14, 14, 13, 13)]
        # type1Std = [(0, 12, 12, 11, 11), (1, 39, 12, 11, 11), (2, 18, 6, 5, 5), (2, 20, 10, 1, 1), (2, 27, 36, 6, 6)]
        type1Std = [(2, 0, 12, 12, 11, 11), (2, 3, 18, 6, 5, 5), (2, 3, 20, 10, 1, 1), (2, 3, 27, 36, 6, 6), (2, 1, 39, 12, 11, 11)]
        # print 'X'

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
            # print 'ret type2'
            return "type2"

    if len(struct) == 3:
        # type1Std = [(0, 24, 20, 9, 9), (1, 10, 10, 9, 9), (2, 15, 19, 4, 4), (2, 30, 32, 5, 5), (3, 14, 14, 13, 13)]
        # type1Std = [(0, 12, 12, 11, 11), (1, 39, 12, 11, 11), (2, 27, 36, 9, 11)]
        type1Std = [(2, 0, 12, 12, 11, 11), (2, 3, 27, 36, 9, 11), (2, 1, 39, 12, 11, 11)]

        return 'type3'
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
            # print 'ret type3'
            return "type3"


    # print 'ret None'
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
    return {k:[(2, 0, 12, 12, 11, 11), (2, 3, 19, 7, 6, 6), (2, 3, 27, 36, 6, 6), (2, 1, 39, 12, 11, 11)] for k in kk}
