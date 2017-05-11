import sys

def toFloat(x):
    try:
        y = float(x)
    except:
        y = -1
    return y

def read(fname):
    reading = False

    d = {}
    with open(fname) as f:
        for line in f:
            line = line.strip()

            if len(line) > 0 and line[0] == '.':
                reading = True

            # because start reading AFTER the line with the '.'
            if reading and line[0] != '.':
                data = line.split('\t')

                # only care about percentages, not ranks
                d[str(data[0])] = map(toFloat, data[1::2])
    return d

def findDiff(d1, d2):
    diffD = {}

    maxDiff = 0
    for shape, d1Vals in d1.iteritems():
        d2Vals = d2[shape]

        # print shape,':'
        # print d1Vals
        # print d2Vals


        # differences:
        diff = [abs(x-y) for (x,y) in zip(d1Vals, d2Vals)]
    
        # print "shape =", shape
        # print diff
        diffD[shape] = diff
        maxDiff = max(maxDiff, max(diff))

    return diffD, maxDiff


def main():
    # arg1 = tbl file 1
    # arg2 = tbl file 2
    if len(sys.argv) < 3:
        sys.stderr.write("Need 2 filenames as arguments!\n")
        return

    f1name = sys.argv[1]
    f2name = sys.argv[2]

    d1 = read(f1name)
    d2 = read(f2name)

    diff, maxDiff = findDiff(d1,d2)

    # for k,v in diff.iteritems():
    #     print k, ':'

    #     # print d1[k]
    #     # print d2[k]
    #     print v
    print maxDiff








if __name__ == "__main__":
    main()