def structTypeYeast(struct):
    introns = ""
    helix = ""
    level0 = [(l1,l2,i,j,w1,w2) for (l1,l2,i,j,w1,w2) in struct if (l1,l2) == (1,0)]
    level1 = [(l1,l2,i,j,w1,w2) for (l1,l2,i,j,w1,w2) in struct if (l1,l2) == (1,2)]
    level2 = [(l1,l2,i,j,w1,w2) for (l1,l2,i,j,w1,w2) in struct if (l1,l2) == (3,2)]

    if level0 == [] and level2 == []:
        introns = "I1, I2 detached"
    elif  level0 == [] or level0 == [(1, 0, 6, 10, 1, 1)]:
        introns = "I1 detached"
    elif  level2 == []:
        introns = "I2 detached"
    else:
        introns = "I1, I2 attached"

    has1A = False
    has1B = False
    has1A_B = False
    for win in level1:
        if win == (1, 2, 16, 18, 4, 4):
            has1A = True
        if win == (1, 2, 19, 23, 2, 2):
            has1B = True 
        if win == (1, 2, 19, 23, 7, 9):
            has1A_B = True


    if has1B and has1A:
        helix = "both helices"
    elif has1B:
        helix = "1B"
    elif has1A:
        helix = "1A"
    elif has1A_B:
        helix = "1AB_joint"
    else:
        helix = "no helix"

    return (introns, helix)


def legend():

    ll = ['I1 detached', 'I1, I2 attached']
    jj = ['both helices', '1A', '1B', '1AB_joint']
    kk = list(itertools.product(ll, jj))
    # kk = [('I1 detached', '1A') ,('I1, I2 attached', '1A') ,('I1, I2 attached', 'both helices'),('I1 detached', 'both helices'),('I1 detached', '1AB_joint') ,('I1, I2 attached', '1AB_joint')]
    legend_ = {a:a for a in kk}

    return legend_

