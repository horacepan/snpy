import os
import random
import pickle
import sys
import pdb
import math
import itertools
import time
from utils import check_memory, partitions
import numpy as np
from young_tableau import YoungTableau, FerrersDiagram
from perm2 import sn

# TODO: make this a tiered dict?
YOR_CACHE = {}
YOR_T_CACHE = {}
CACHE = {'hit': 0, 'sparse_hit': 0}

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
    if (ferrers, permutation.tup_rep) in YOR_CACHE:
        CACHE['hit'] += 1
        return YOR_CACHE[(ferrers.partition, permutation.tup_rep)]

    if all(map(lambda x: len(x) <= 1, permutation.cycle_decomposition)):
        # TODO: make a static/class function for this
        n = len(FerrersDiagram.TABLEAUX_CACHE[ferrers.partition])
        YOR_CACHE[(ferrers.partition, permutation.tup_rep)] = np.eye(n)
        return YOR_CACHE[(ferrers.partition, permutation.tup_rep)]

    res = None
    for cycle in permutation.cycle_decomposition:
        ts = []
        for t in cycle_to_adj_transpositions(cycle, ferrers.size):
            y = yor_trans(ferrers, t)
            ts.append(t)
            if res is None:
                res = y
            else:
                res = res.dot(y)
    if use_cache:
        YOR_CACHE[(ferrers.partition, permutation.tup_rep)] = res
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
        CACHE['hit'] += 1
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

def ysemi(ferrers, permutation):
    '''
    Compute the irreps of the given shape using Young's Seminormal Form

    ferrers: FerrersDiagram
    partition: tuple of ints
    permutation: list of tuples in standard cycle notation. Order doesnt matter
                 since the cycles are disjoint in standard notation
    Returns: an irrep matrix of size d x d, where d is the number of standard tableaux of the
    given FerrersDiagram shape
    '''
    # TODO: This is a hacky way of passing in the identity permutation
    if all(map(lambda x: len(x) <= 1, permutation.cycle_decomposition)):
        return np.eye(len(FerrersDiagram.TABLEAUX_CACHE[ferrers.partition]) )

    res = None
    for cycle in permutation.cycle_decomposition:
        # rewrite the cycle in terms of the adjacent transpositions group generators
        for t in reversed(cycle_to_adj_transpositions(cycle, ferrers.size)):
            if t[0] > t[1]:
                # keep transpositions in standard notation
                t = (t[1], t[0])
            y = ysemi_t(ferrers, t)
            if res is None:
                res = y
            else:
                res = y.dot(res)

    return res

def ysemi_t(f, transposition):
    '''
    Young's seminormal form for adjacent transposition permutations.
    EG: permutations of the form (k, k+1)

    f: ferrers diagram
    Returns a matrix of size = n_p x n_p, where
    n_p = the number of young tableau of the given partition
    '''
    tabs = f.tableaux
    rep = np.zeros((len(tabs), len(tabs)))

    for t in tabs:
        res = t.transpose(transposition)
        if res is None:
            # see if same row/col
            # same row --> 1
            if t.get_row(transposition[0]) == t.get_row(transposition[1]):
                rep[t.idx, t.idx] = 1.
            # same col --> -1
            elif t.get_col(transposition[0]) == t.get_col(transposition[1]):
                rep[t.idx, t.idx] = -1.
            continue

        i = t.idx
        j = res.idx
        if i < j:
            dist = t.ax_dist(*transposition)
            rep[i, i] = 1. / dist
            # rep[i, j] = 1. - (1./(dist**2))
            rep[i, j] = 1. + (1.0 / dist)
            rep[j, i] = 1.
            rep[j, j] = -1. / dist
            rep[j, i] = 1. - (1.0 / dist)

    return rep

# TODO: Benchmarking function should go elsewhere
def benchmark(n):
    '''
    Benchmark time/memory usage for generating all YoungTableau for S_8
    '''
    tstart = time.time()
    _partitions = partitions(n)
    s_n = sn(n)
    print('Starting...')
    times = []
    for idx, p in enumerate(_partitions):
        start = time.time()
        sdict = {}
        f = FerrersDiagram(p)
        if os.path.exists('/local/hopan/irreps/s_9/{}.pkl'.format(p)):
            print('Skipping {}'.format(p))
            continue
        for perm in s_n:
            start = time.time()
            y = yor(f, perm)
            end = time.time()
            if random.random() > 0.1 and len(times) < 1000:
                times.append(end - start)
            if len(times) >= 1000:
                pdb.set_trace()
            sdict[perm.tup_rep] = y

        done = time.time() - start
        print('Elapsed: {:.2f}mins | Done {} / {} | Partition: {}'.format(done / 60., idx, len(_partitions), p))

        with open('/local/hopan/irreps/s_{}/{}.pkl'.format(n, p), 'wb') as f:
            pickle.dump(sdict, f, protocol=pickle.HIGHEST_PROTOCOL)

    tend = time.time() - tstart
    print('Total time compute yor matrices for S_{}: {:3f}'.format(n, tend))
    print(CACHE)

def load_yor(fname, partition):
    #print('loading yor from: {}'.format(fname))
    with open(fname, 'rb') as f:
        yor_dict = pickle.load(f)
        # mapping form permutation in list form to numpy array
        for perm, mat in yor_dict.items():
            YOR_CACHE[(partition, perm)] = mat

        return yor_dict

if __name__ == '__main__':
    n = int(sys.argv[1])
    benchmark(n)
