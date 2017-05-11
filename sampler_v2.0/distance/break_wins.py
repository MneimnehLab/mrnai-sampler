#!/usr/bin/python

'''
tool to break sampler output windows into smaller windows of equal size based on sequence alignment

a new window starts whenever there is a gap on one level... so with bonds (1,1) and (1,2), 
we have two windows of size 1
'''

from collections import defaultdict
import sys

class WinBreaker(object):

    def __init__(self, seqs, S):
        self.seqs = seqs
        self.S = S
        self.allowed_pairs = [('A','U'), ('U','A'), ('C','G'), ('G','C'), ('G','U'), ('U','G')] 


    def align_bases(self, win):
        (l1,l2, i, j, w1, w2) = win
        
        s1 = self.seqs[l1]
        s2 = self.seqs[l2][::-1] # odd sequence is reversed
    
        seq1 = "." + s1[i-w1-1:i]
        seq2 = "." + s2[j-w2-1:j]

        # print "finding LCS of ", seq1, seq2

        
        m, n = len(seq1), len(seq2)
        # print 'm,n =', m,n
        H = defaultdict(int)
        backtrack = {}
        for i in range(1,m):
            for j in range(1,n):
                # print i,j, seq1[i],seq2[j]
                # if seq1[i] == seq2[j]:
                if (seq1[i], seq2[j]) in self.allowed_pairs:
                    H[(i,j)] = H[(i-1,j-1)] + 1
                    backtrack[(i,j)] = (i-1,j-1)
                elif H[(i,j-1)] > H[(i-1,j)]:
                    backtrack[(i,j)] = (i,j-1)
                    H[(i,j)] = H[(i,j-1)]
                else:
                    backtrack[(i,j)] = (i-1,j)
                    H[(i,j)] = H[(i-1,j)]
                

        # print H[(m-1,n-1)]
        i = m-1
        j = n-1
        bonds = []
        while i > 0 and j > 0:
            t = backtrack[(i,j)]
            # print 't =', t
            if i==t[0]+1 and j==t[1]+1:
                # print i,j
                # re-add original indices
                bonds.append( (win[2]-win[4]-1+i, win[3]-win[5]-1+j )  )
            i = t[0]
            j = t[1]

        return bonds[::-1]


    def break_wins(self):
        new_wins = []
        
        
        for win in self.S:
            prev_i, prev_j = -100, -100
            bonds = self.align_bases(win)
            # print bonds
            # new window starts if a gap is detected on one level
            bond_list = []
            for (i,j) in bonds:
                if i == prev_i+1 and j == prev_j+1:
                    # same window
                    bond_list.append((i,j))
                else:
                    if bond_list != []:
                        all_i = sorted([i1 for (i1,j1) in bond_list])
                        all_j = sorted([j1 for (i1,j1) in bond_list])
                        window = (win[0], win[1], all_i[-1], all_j[-1], all_i[-1] - all_i[0], all_j[-1] - all_j[0]) 
                        new_wins.append(window)
            
                    bond_list = [(i,j)]
                prev_i = i
                prev_j = j

            if bond_list != []:
                all_i = sorted([i1 for (i1,j1) in bond_list])
                all_j = sorted([j1 for (i1,j1) in bond_list])
                window = (win[0], win[1], all_i[-1], all_j[-1], all_i[-1] - all_i[0], all_j[-1] - all_j[0]) 
                new_wins.append(window)

        return new_wins



def main():
    
    if len(sys.argv) == 1:
        sys.stderr.write('Need sequences file as argument!\n\n')
        return

    try:
        with open(sys.argv[1]) as f:
            seqs = f.readline().strip().split('&')
            if len(seqs) < 2: # should contain at least two sequences
                raise
            # for s in seqs:
            #     sys.stderr.write('Read sequence ' + s + '\n')
    except:
        sys.stderr.write('Sequences file doesn\'t exist or is malformed!\n\n')
        return


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
        except:
            sys.stderr.write('Error parsing line ' + str(line_num+1) + '!\n\n')
            return

        try:
            wb = WinBreaker(seqs, wlist)
            new_struct =  wb.break_wins()

            if len(sp) >= 3:
                output = sp[0] + '\t' + str(new_struct) + '\t' + '\t'.join(sp[2:])
            else:
                output = sp[0] + '\t' + str(new_struct)
            
            print output

        except Exception as e:
            sys.stderr.write('Error in parsing... sequences data is probably invalid... fewer or shorter sequences\n')
            sys.stderr.write('... program will terminate.\n\n')
            # sys.stderr.write('Actual error: \n')
            # sys.stderr.write('\t' + str(type(e)) + '\n')
            # sys.stderr.write('\t' + str(e.args) + '\n')
            # sys.stderr.write('\t' + str(e) + '\n')
            return


if __name__ == '__main__':
    main()






def test():
    seq1 = 'CGGUUUAAGUGGGCCCCGGUAAUCUUUUCGUACUCGCCAAAGUUGAAGAAGAUUAUCGGGGUUUUUGCUU'
    seq2 = 'AAGCAAAAACCCCGAUAAUCUUCUUCAACUUUGGCGAGUACGAAAAGAUUACCGGGGCCCACUUAAACCG'


    # unit tests
    test_structs = [
        [(0,1, 12, 13, 10, 11), (0,1, 50, 50, 23, 24), (0,1, 69, 69, 3, 3)],
        [(0,1, 15, 15, 12, 12), (0,1, 50, 50, 24, 24), (0,1, 69, 69, 7, 7)],
        [(0,1, 1,1,0,0), (0,1,2,3,0,0), (0,1,5,5,1,1)]    
    ]

    for S in test_structs:
        wb = WinBreaker([seq1, seq2], S)
        print wb.break_wins()

