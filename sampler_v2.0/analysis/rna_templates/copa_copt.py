sequences = 'CGGUUUAAGUGGGCCCCGGUAAUCUUUUCGUACUCGCCAAAGUUGAAGAAGAUUAUCGGGGUUUUUGCUU&AAGCAAAAACCCCGAUAAUCUUCUUCAACUUUGGCGAGUACGAAAAGAUUACCGGGGCCCACUUAAACCG'.split('&')

def legend():
    return {"A":"first, middle", "B":"first split, middle", "C":"first,  middle split", "D":"first, middle, last", "E":"first split, middle, last", "F":"first,  middle split, last"}

def rfrncForClass():
    return {k:[(0, 1, 12, 12, 11, 11), (0, 1, 35, 35, 8, 8), (0, 1, 47, 47, 5, 5)] for k in ['A','B','C','D','E','F']}
        

def struct_type(struct):

    MIDDLE_RIGHT_WIN_END = 50
    MIDDLE_LEFT_WIN_START = 24

    if len(struct) == 2:
        A = struct[0]
        B = struct[1]
        
        # first, middle
        if                                         A[2] <= 15 and A[3] <= 15 and\
           B[2]-B[4] >= 24 and B[3]-B[5] >= 24 and B[2] <= 50 and B[4] <= 50:
            return "A"

    elif len(struct) == 3:
        A = struct[0]
        B = struct[1]
        C = struct[2]

        # first split, middle
        if                                         B[1] <= 15 and B[2] <= 15 and    \
           C[2]-C[4] >= 24 and C[3]-C[5] >= 24 and C[2] <= 50 and C[3] <= 50:
            return "B"      
          
        
        #first,   middle split
        elif A[2] <= 15 and A[3] <= 15 and \
            B[2]-B[4] >= 25 and B[3]-B[5] >= 25 and B[2] <= MIDDLE_RIGHT_WIN_END and B[3] <= MIDDLE_RIGHT_WIN_END and\
            C[2]-C[4] >= MIDDLE_LEFT_WIN_START and C[3]-C[5] >= MIDDLE_LEFT_WIN_START and C[2] <= 50 and C[3] <= 50:
            return "C"

            
        # first, middle, last
        elif                                         A[2] <= 15 and A[3] <= 15 and\
             B[2]-B[4] >= 24 and B[3]-B[5] >= 24 and B[2] <= 50 and B[3] <= 50 and\
             C[2]-C[4] >= 60 and C[3]-C[5] >= 60 and C[2] <= 70 and C[3] <= 70:
            return "D"

    elif len(struct) == 4:
        A = struct[0]
        B = struct[1]
        C = struct[2]
        D = struct[3]

        # first split, middle, last
        if                                         B[1] <= 15 and B[2] <= 15 and \
           C[2]-C[4] >= 24 and C[3]-C[5] >= 24 and C[2] <= 50 and C[3] <= 50 and \
           D[2]-D[4] >= 60 and D[3]-D[5] >= 60 and D[2] <= 70 and D[3] <= 70:
            return "E"

        # first,  middle split, last
        if                                         A[1] <= 15 and A[2] <= 15 and \
           B[2]-B[4] >= 24 and B[3]-B[5] >= 24 and B[2] <= MIDDLE_RIGHT_WIN_END and B[3] <= MIDDLE_RIGHT_WIN_END and \
           C[2]-C[4] >= 39 and C[3]-C[5] >= 39 and C[2] <= 50 and C[3] <= 50 and \
           D[2]-D[4] >= 60 and D[3]-D[5] >= 60 and D[2] <= 70 and D[3] <= 70:
            
            return "F"

    return None


def legend():
    return {"A":"first, middle", "B":"first split, middle", "C":"first,  middle split", "D":"first, middle, last", "E":"first split, middle, last", "F":"first,  middle split, last"}


if __name__ == '__main__':
    S = [(0, 1, 14, 14, 11, 11), (0, 1, 38, 38, 11, 11), (0, 1, 48, 48, 7, 8)]

    print struct_type(S)

