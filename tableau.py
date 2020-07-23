from utils import *
from young_tableau import FerrersDiagram, get_minus_partition
import time
import pdb
import copy
from functools import total_ordering

def swap(x, i, j):
    if x == i:
        return j
    elif x == j:
        return i
    else:
        return x

@total_ordering
class YoungTableau:
    CACHE = {}
    def __init__(self, contents, partition, n):
        self.size = n
        self.partition = partition
        self.contents = contents
        self._val = tuple(i for row in contents for i in row)
        self._rows = {}
        self._cols = {}
        self._fill_rows_cols()
        self._idx = None

        if partition not in YoungTableau.CACHE:
            YoungTableau.CACHE[partition] = {}
        YoungTableau.CACHE[partition][self._val] = self

    @property
    def index(self):
        return self._idx

    def set_index(self, idx):
        self._idx = idx

    def get_row(self, x):
        return self._rows[x]

    def get_col(self, x):
        return self._cols[x]

    def ax_dist(self, i, j):
        ri = self.get_row(i)
        rj = self.get_row(j)
        ci = self.get_col(i)
        cj = self.get_col(j)

        return (cj - ci) - (rj - ri)

    def transpose(self, a, b):
        '''
        Returns the Young Tableau that results when you apply the given transposition (a, b)
        to this tableau.
        '''
        ra = self.get_row(a)
        rb = self.get_row(b)
        ca = self.get_col(a)
        cb = self.get_col(b)

        if ra != rb and ca != cb:
            swapped_tab = tuple(swap(x, a, b) for x in self._val)
            return YoungTableau.CACHE[self.partition].get(swapped_tab, None)
        else:
            # (a, b) * this != valid
            return None

    def _fill_rows_cols(self):
        for ridx, row in enumerate(self.contents):
            for cidx, val in enumerate(row):
                self._rows[val] = ridx
                self._cols[val] = cidx

    def __repr__(self):
        rep_str = ''
        for idx, l in enumerate(self.contents):
            rep_str += ''.join(map(lambda x: '[{}]'.format(x), l))
            if idx != len(self.contents) - 1:
                rep_str += '\n'
        return rep_str

    def __lt__(self, other):
        for i in range(self.size, 0, -1):
            # get row of i in self and in other
            if self.get_row(i) < other.get_row(i):
                return True
            elif self.get_row(i) > other.get_row(i):
                return False
            # otherwise they're in the same row and you proceed
        return False

def gen_standard_tableaux(partition):
    n = sum(partition)
    tabs = [[]]

    for i in range(n):
        curr_tabs = []
        for t in tabs:
            nts = _gen_next(t, partition, i + 1)
            curr_tabs.extend(nts)

        tabs = curr_tabs

    return tabs

def _gen_next(curr, partition, el):
    next_parts = []

    for idx, row in enumerate(curr):
        if len(row) < partition[idx]: # there is space on this row
            if idx == 0 or (idx > 0 and len(row) < len(curr[idx - 1])):
                new_tab = copy.deepcopy(curr)
                new_tab[idx].append(el)
                next_parts.append(new_tab)

    if len(curr) < len(partition):
        new_tab = copy.deepcopy(curr)
        new_tab.append([el])
        next_parts.append(new_tab)

    return next_parts[::-1]

# TODO: Should just generate the standard tableaux in sorted order to avoid post sorting
def compute_syt(partition):
    n = sum(partition)
    yts = sorted([YoungTableau(t, partition, n) for t in gen_standard_tableaux(partition)])
    for idx, tab in enumerate(yts):
        tab.set_index(idx)

    return yts

def pp(tableaux):
    for row in tableaux:
        fmt = '[{}]' * len(row)
        rowstr = fmt.format(*row)
        print(rowstr)

def check_eq(yt, ft):
    for ry, rf in zip(yt.contents, ft.contents):
        if tuple(ry) != rf:
            return False
    return True

def branch_down(partition):
    '''
    Return a list of FerrersDiagrams of all the partitions
    you'd get if you removed a corner node
    '''
    if sum(partition) == 1:
        return []

    branches = []
    for idx in range(len(partition)):
        # can branch if partition[idx] - 1 >= partition[idx + 1]
        next_row_len = 0
        if idx != len(partition) - 1:
            next_row_len = partition[idx + 1]
        if partition[idx] - 1 >= next_row_len:
            # valid
            new_partition = get_minus_partition(partition, idx)
            branches.append(new_partition)

    return branches


if __name__== '__main__':
    st = time.time()
    part = (28, 1, 1)

    n = sum(part)
    sty = gen_standard_tableaux(part)
    sty_time = time.time() - st

    st = time.time()
    yts = [YoungTableau(t, part, n) for t in sty]

    yt_time = time.time() - st
    tot = sty_time + yt_time

    st = time.time()
    yts = sorted(yts)
    end = time.time()
    sort_time = end - st
    tot += sort_time
    print('Num sty', len(yts))
    print(f'Elapsed: {tot:.2f} | Gen Shapes: {sty_time:.2f}s | Gen YT Obj: {yt_time:.2f}s')
    phaperint('Sort time: {:.2f}s'.format(sort_time))
