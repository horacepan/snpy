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
