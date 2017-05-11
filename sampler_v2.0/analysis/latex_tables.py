import sys


def main():

    percents = {}
    ranks = {}
    f1s = {}

    tbls = {}
    tbls['percent'] = {}
    tbls['ranks'] = {}
    tbls['f1'] = {}

    distinctLabels = []  # not using set() because we want to maintain original order 

    reading = None
    for line in sys.stdin:
        line = line.strip()
        if len(line) < 3:
            continue


        if line[:3] == 'per':
            reading = 'percent'
            continue
        elif line[:3] == 'avg':
            reading = 'ranks'
            continue
        elif line[:3] == 'f1s':
            reading = 'f1'
            continue

        # print 'for', reading
        # print line
        data = line.split()
        # print data
        label = data[0]
        vals = map(float, data[1:])
        # print label
        # print vals
        
        if label not in distinctLabels:
            distinctLabels.append(label)

        # sanitize vals to replace 0's with '-'s
        for i in range(len(vals)):
            if vals[i] == 0:
                if reading != 'percent':
                    vals[i] = '-'
        tbls[reading][label] = vals


    # for k,v in tbls.iteritems():
    #     print "for", k

    #     for k2,v2 in v.iteritems():
    #         print k2,v2


    legend = {'type1':'big window', 'type2':'equal split', 'type3':'unequal split'}

    
    print r'\begin{table*}'
    print r'\begin{center} '
    print r'\begin{tabular}{l|ccc|ccc|ccc|ccc|ccc} '
    print r'k & 1 &&& 2 &&& 3 &&& 4 &&& 5\\ \hline'

    
    for label in distinctLabels:
        row = []
        # print label, '\t',
        row.append(legend[label])
        for num in range(5):
            # print percetnages, rank, f1
            # print tbls['percent'][label][num],'\t',
            # print tbls['ranks'][label][num],'\t',
            # print tbls['f1'][label][num],'\t',
            row.append(str(tbls['percent'][label][num]))
            row.append(str(tbls['ranks'][label][num]))
            row.append(str(tbls['f1'][label][num]))
        # print ''
        print ' & '.join(row), 

        if label == distinctLabels[-1]:
            print r'\vspace{1em}',
        print r'\\'

    print r'k & 6 &&& 7 &&& 8 &&& 9 &&& 10\\ \hline'

    for label in distinctLabels:
        row = []
        # print label, '\t',
        row.append(legend[label])
        for num in range(5,10):
            # print percetnages, rank, f1
            # print tbls['percent'][label][num],'\t',
            # print tbls['ranks'][label][num],'\t',
            # print tbls['f1'][label][num],'\t',
            row.append(str(tbls['percent'][label][num]))
            row.append(str(tbls['ranks'][label][num]))
            row.append(str(tbls['f1'][label][num]))
        # print ''
        print ' & '.join(row), 
        print r'\\'


    print r'\vspace{1em}'
    print r'\end{tabular}'
    print r'\caption{Caption} '
    print r'\end{center}'
    print r'\end{table*}'

    print ''
    print ''
    print ''
    print ''
    print ''
    print ''
    print ''



if __name__ == '__main__':
    main()

