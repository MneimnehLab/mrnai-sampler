def struct_type(struct):

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


def legend():
    return {"A":"first, middle", "B":"first split, middle", "C":"first,  middle split", "D":"first, middle, last", "E":"first split, middle, last", "F":"first,  middle split, last"}
