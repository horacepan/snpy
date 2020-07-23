import numpy as np
from young_tableau import FerrersDiagram
from perm import sn

YOR_CACHE = {}
YOR_T_CACHE = {}

def cycle_to_adj_transpositions(cyc, n):
    '''
    cyc: tuple of ints, the permutation cycle
    n: size of the S_n group that this permutation lives in
    Given a cycle: a -> b -> c
    Generate perm = [cyc(i) for i in range(1, n+1)]
    TODO: can we do this without creating the mapping list
    '''
    # do bubble sort to turn this into a product of adjacent transpositions
    cyc_map = lambda x: x if x not in cyc else cyc[(cyc.index(x) + 1) % len(cyc)]
    perm = [ cyc_map(i) for i in range(1, n+1)]
    factors = []

    for i in range(n):
        for j in range(n-1, i, -1):
            if perm[j] < perm[j-1]:
                perm[j], perm[j-1] = perm[j-1], perm[j]
                factors.append((j, j+1))

    return list(reversed(factors))

def perm_to_adj_transpositions(perm, n):
    '''
    perm: a permutation of S_n list of tuple of ints. perm is given in its canonical cycle decomposition
    format
    n: integer, size of the permutation group that perm is a member of
    '''
    all_trans = []
    for cyc in perm:
        all_trans.extend(cycle_to_adj_transpositions(cyc, n))

    return all_trans

def yor(ferrers, permutation, use_cache=True):
    '''
    Compute the irreps of a given shape using Young's Orthogonal Representation (YOR)

    ferrers: FerrersDiagram
    permutation: perm.Perm object
    permutation: list of tuples for the permutation in disjoint
                 cycle notation
    Returns: an irrep matrix of size d x d, where d is the number of standard tableaux of the
    given FerrersDiagram shape
    '''
    if (ferrers, permutation.tup) in YOR_CACHE:
        return YOR_CACHE[(ferrers.partition, permutation.tup)]

    if all(map(lambda x: len(x) <= 1, permutation.cycle_decomposition)):
        # TODO: make a static/class function for this
        n = len(FerrersDiagram.TABLEAUX_CACHE[ferrers.partition])
        YOR_CACHE[(ferrers.partition, permutation.tup)] = np.eye(n)
        return YOR_CACHE[(ferrers.partition, permutation.tup)]

    res = None
    all_ts = []
    for cycle in permutation.cycle_decomposition:
        ts = []
        for t in cycle_to_adj_transpositions(cycle, ferrers.size):
            y = yor_trans(ferrers, t)
            ts.append(t)
            if res is None:
                res = y
            else:
                res = res.dot(y)
        all_ts.append(ts)

    if use_cache:
        YOR_CACHE[(ferrers.partition, permutation.tup)] = res
    return res

def yor_trans(ferrers, transposition):
    '''
    Young's seminormal form for adjacent transposition permutations.
    EG: permutations of the form (k, k+1)
    ferrers: a FerrersDiagram object
    transposition: a 2-tuple of ints

    Returns: an irrep matrix of size d x d, where d is the number of standard tableaux of the
    given FerrersDiagram shape
    '''
    assert transposition[0] < transposition[1]
    if (ferrers, transposition) in YOR_T_CACHE:
        return YOR_T_CACHE[(ferrers, transposition)]

    tabs = ferrers.tableaux
    rep = np.zeros((len(tabs), len(tabs)))
    for tab in tabs:
        other = tab.transpose(transposition)
        dist = tab.dist(*transposition)
        i = tab.idx
        rep[i, i] = 1. / dist

        if other is not None:
            j = other.idx
            rep[i, j] = np.sqrt(1 - (1. / dist) ** 2)
            rep[j, i] = rep[i, j]
            rep[j, j] = 1. / other.dist(*transposition)

    YOR_T_CACHE[(ferrers, transposition)] = rep
    return rep
