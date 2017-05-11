#!/usr/bin/python

from collections import defaultdict
import sys

INF = 999999

# winsorter = lambda x,y : (x[0],x[1]) > (y[0],y[1])

def get_boundary(S, m):
    pegs_by_level = defaultdict(list)
    max_level = 0

    for (l1,l2,i,j,w1,w2) in S:
        pegs_by_level[l1].extend(range(i-w1,i+1))
        pegs_by_level[l2].extend(range(j-w2,j+1))
        max_level = max(max_level, l1, l2)

    boundary = [-INF] * m

    for level, pegs in pegs_by_level.iteritems():
        if len(pegs) > 0:
            boundary[level] = max(pegs)

    return boundary


def get_terminal_win(S, m):
    if len(S) == 0:
        return None

    boundary = get_boundary(S, m)

    on_boundary = []

    for (l1,l2,i,j,w1,w2) in S:
        if boundary[l1] == i and boundary[l2] == j:
            on_boundary.append((l1,l2,i,j,w1,w2))

    on_boundary.sort()

    # print 'on_boundary =', on_boundary
    return on_boundary[-1]


def sort_structure(S):
    new_set = []

    m = 0
    for (l1,l2,i,j,w1,w2) in S:
        m = max(m, l1,l2)

    m += 1

    while S:
        term = get_terminal_win(S, m)
        new_set.append(term)
        S = [w for w in S if w != term]

    return new_set[::-1]




def main():
    
    for line_num,line in enumerate(sys.stdin):        
        try:
            sp = line[:-1].split("\t")
            
            # handle both file formats
            if len(sp) >= 3:
                num, wlist, energy = sp[:3]
            elif len(sp) == 1:
                wlist = line.strip()
            else:
                raise
            
            wlist = eval(wlist)
            if type(wlist) is not list:
                raise    

            sorted_struct = sort_structure(wlist)
            if len(sp) >= 3:
                output = sp[0] + '\t' + str(sorted_struct) + '\t' + '\t'.join(sp[2:])
            else:
                output = sp[0] + '\t' + str(sorted_struct)
            
            print output

        except Exception as e:
            sys.stderr.write('Error sorting line ' + str(line_num+1) + '!\n\n')
            return

        




if __name__ == "__main__":
    main()



