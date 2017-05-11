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
                vals[i] = '-'
        tbls[reading][label] = vals


    # for k,v in tbls.iteritems():
    #     print "for", k

    #     for k2,v2 in v.iteritems():
    #         print k2,v2



    for label in distinctLabels:
        print label, '\t',
        for num in range(10):
            # print percetnages, rank, f1
            print tbls['percent'][label][num],'\t',
            print tbls['ranks'][label][num],'\t',
            print tbls['f1'][label][num],'\t',
        print ''




if __name__ == '__main__':
    main()

