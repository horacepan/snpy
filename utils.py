import math

def partitions(n, start=1):
    '''
    Generate all the partitions of n
    n: integer
    start: integer
    Returns a list of int tuples
    '''
    if n == 0:
        return [()]

    parts = [(n, )]
    for i in range(start, n // 2 + 1):
        for p in partitions(n-i, i):
            parts.append(p + (i, ))

    return parts

def hook_length(partition):
    n = sum(partition)
    numerator = math.factorial(n)
    denominator = 1

    for idx, rowlen in enumerate(partition):
        for _idx in range(rowlen):
            rl = rowlen - _idx
            cl = len([r for r in partition[idx+1:] if r > _idx])
            denominator *= (rl + cl)

    return numerator / denominator

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

def cycle_decomposition(perm_tup):
    pass
