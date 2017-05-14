sequences = 'AAAGAGAGAAGUGAACCAGAGAAACACACGCG&UCACAGUCCUCUUU&CGCGUGCUAUAUUACCUGGUA'.split('&')

def legend():
    # kk = ["type1", "type2"]
    kk = ["type"+str(t) for t in range(1,3+1)]
    # kk = [('I1 detached', '1A') ,('I1, I2 attached', '1A') ,('I1, I2 attached', 'both helices'),('I1 detached', 'both helices'),('I1 detached', '1AB_joint') ,('I1, I2 attached', '1AB_joint')]
    return {a:a for a in kk}

def rfrncForClass():
    kk = ["type"+str(t) for t in range(1,3+1)]
    # 
    # return {k:[(0, 6, 6, 5, 5), (0, 14, 14, 3, 3),  (1, 19, 6, 4, 4), (1, 32, 21, 5, 5)] for k in kk}
    return {k:[(0, 1, 6, 6, 5, 5), (0, 1, 14, 14, 3, 3), (0, 2, 19, 6, 4, 4), (0, 2, 32, 21, 5, 5)] for k in kk}


def struct_type(struct):
    # print 's =', struct
    DELTA = 3

    if len(struct) == 3:
        # type1Std = [(0, 24, 20, 9, 9), (1, 10, 10, 9, 9), (2, 15, 19, 4, 4), (2, 30, 32, 5, 5), (3, 14, 14, 13, 13)]
        # type1Std = [(0, 14, 14, 13, 13),  (1, 19, 6, 4, 4), (1, 32, 21, 5, 5)]
        type1Std = [(0, 1, 14, 14, 13, 13), (0, 2, 19, 6, 4, 4), (0, 2, 32, 21, 5, 5)]


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
            return "type1"

    elif len(struct) == 4:
        # type2Std = [(0, 24, 20, 9, 9), (1, 10, 10, 9, 9), (2, 15, 19, 4, 4), (2, 30, 32, 5, 5), (3, 6, 6, 5, 5), (3, 14, 14, 3, 3)]
        # type2Std =  [(0, 6, 6, 5, 5), (0, 14, 14, 3, 3),  (1, 19, 6, 4, 4), (1, 32, 21, 5, 5)]
                    # [(0, 9, 8, 8, 7), (0, 14, 14, 3, 3), (1, 19, 6, 4, 4), (1, 27, 18, 2, 2), (1, 32, 21, 2, 2)]
        type2Std = [(0, 1, 6, 6, 5, 5), (0, 1, 14, 14, 3, 3), (0, 2, 19, 6, 4, 4), (0, 2, 32, 21, 5, 5)]

        # print 's4 =', struct
        
        
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
            if struct[0][4] == struct[0][5]:
                return "type2"
            else:
                # print struct
                return "type3"

    elif len(struct) == 5:
        # type2Std =  [(0, 9, 8, 8, 7), (0, 14, 14, 3, 3), (1, 19, 6, 4, 4), (1, 27, 18, 2, 2), (1, 32, 21, 2, 2)]
        type2Std = [(0, 1, 6, 6, 5, 5), (0, 1, 14, 14, 3, 3), (0, 2, 19, 6, 4, 4), (0, 2, 27, 18, 2, 2), (0, 2, 32, 21, 2, 2)]
        
        
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
            if struct[0][4] != struct[0][5]:
                # print 'X'
                return "type3"
            # else:
            #     # print struct
            #     return "type3"




    return None


def main():
    # print struct_type([(0, 1, 9, 8, 8, 7), (0, 1, 14, 14, 3, 3), (0, 2, 19, 6, 4, 4), (0, 2, 32, 21, 5, 5)])
    pass


if __name__ == '__main__':
    main()

