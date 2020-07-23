import os
import resource
import pdb
import time
import itertools
import numpy as np
from functools import total_ordering

FERRERS_CACHE = {}
def swap(x, i, j):
    if x == i:
        return j
    elif x == j:
        return i
    else:
        return x

def set_idx(sorted_tabs):
    for idx, tab in sorted_tabs:
        tab.set_idx(idx)

def get_minus_partition(partition, idx):
    '''
    Return the partition that would result if you subtract 1 from
    the given index.
    Precondition: the input is valid
    Ex:
        (4, 2, 1), 0 --> (3, 2, 1)
        (4, 1, 1), 1 --> (4, 1)
        (4, 2, 1), 2 --> (4, 2)
    '''
    res = []
    for _idx, i in enumerate(partition):
        if idx == _idx and i > 1:
            res.append(i - 1)
        elif idx != _idx:
            res.append(i)
    return tuple(res)

class FerrersDiagram:
    TABLEAUX_CACHE = {}
    def __init__(self, partition):
        '''
        partition: tuple of ints
        '''
        self.partition = partition
        self.size = sum(partition)

        if partition not in FerrersDiagram.TABLEAUX_CACHE:
            self.tableaux = self.gen_tableaux() # list of sorted young tableaux
        else:
            self.tableaux = FerrersDiagram.TABLEAUX_CACHE[partition]
        FERRERS_CACHE[partition] = self

    @staticmethod
    def from_partition(partition):
        if partition not in FERRERS_CACHE:
            FERRERS_CACHE[partition] = FerrersDiagram(partition)
        return FERRERS_CACHE[partition]

    def branch_down(self):
        '''
        Return a list of FerrersDiagrams of all the partitions
        you'd get if you removed a corner node
        '''
        if sum(self.partition) == 1:
            return []

        branches = []
        for idx in range(len(self.partition)):
            # can branch if partition[idx] - 1 >= partition[idx + 1]
            next_row_len = 0
            if idx != len(self.partition) - 1:
                next_row_len = self.partition[idx + 1]
            if self.partition[idx] - 1 >= next_row_len:
                # valid
                new_partition = get_minus_partition(self.partition, idx)
                branches.append(FerrersDiagram.from_partition(new_partition))

        return branches

    def __repr__(self):
        rep_str = ''
        for idx, size in enumerate(self.partition):
            rep_str += '[ ]' *size
            if idx != len(self.partition) - 1:
                rep_str += '\n'
        return rep_str

    def pp_tableaux(self):
        '''Pretty printer for all tableaux'''
        for ta in f.tableaux:
            print(ta)
            print('======')

    # TODO: A LOT of room here for optimization!
    def gen_tableaux(self, perms=None):
        '''
        Returns: sorted list of the YoungTableau of the given shape/partition.
        '''
        tabs = []
        if sum(self.partition) == 0:
            return tabs

        # 1 is always the top left corner so no need to do this
        #for p in itertools.permutations(range(2, self.size+1)):
        if perms is None:
            perms = [(1, ) + p for p in itertools.permutations(range(2, self.size+1))]

        for p in perms:
            try:
                if YoungTableau.valid_static(self.partition, p):
                    yt = make_young_tableau(self.partition, p)
                    tabs.append(yt)
            except:
                pdb.set_trace()
        tabs.sort()

        # set the index which is necessary for computing YOR matrices
        for idx, tab in enumerate(tabs):
            tab.set_idx(idx)

        FerrersDiagram.TABLEAUX_CACHE[self.partition] = tabs
        return tabs

    def n_tabs(self):
        return len(self.tableaux)

def n_tabs(partition):
    fd = FerrersDiagram.from_partition(partition)
    return fd.n_tabs()

def make_young_tableau(shape, vals):
    '''
    shape: tuple of row size of ferrers diagram. Sum of tuple elements is n
    vals: a permutation of {1, 2, ..., n}
    Return a filled in Young Tableau

    Ex: shape = (2, 1), vals = [1, 2, 3]
        Returns: [1][2]
                 [3]

        shape = (2, 2), vals = [1, 3, 2, 4]
        Returns: [1][3]
                 [2][4]
    '''
    contents = []
    i = 0
    for rowsize in shape:
        contents.append(vals[i:i+rowsize])
        i += rowsize
    return YoungTableau(contents, vals)

@total_ordering
class YoungTableau:
    CACHE = {}
    def __init__(self, contents, vals=None):
        '''
        contents: list of tuples of the filled in ferrers blocks
        vals: the full tuple for the permutation
        '''
        self.partition = tuple(map(lambda x: len(x), contents))
        self.contents = contents
        self.size = sum(self.partition)
        self.vals = vals
        self.idx = None

        self._row = {}
        self._col = {}
        for idx, row in enumerate(contents):
            for ri, x in enumerate(row):
                self._row[x] = idx + 1
                self._col[x] = ri + 1

        if self.partition not in YoungTableau.CACHE:
            YoungTableau.CACHE[self.partition] = {}
        YoungTableau.CACHE[self.partition][vals] = self

    def __lt__(self, other):
        '''
        Use last letter ordering to determine which tableau is "larger"
        self is less than other if element n is on a higher row.
        if n is on the same row then check n-1, and so on.
        '''
        for i in range(self.size, 0, -1):
            # get row of i in self and in other
            if self.get_row(i) < other.get_row(i):
                return True
            elif self.get_row(i) > other.get_row(i):
                return False
            # otherwise they're in the same row and you proceed
        return False

    def set_idx(self, idx):
        '''
        Can only know the index of tableaux if you get it
        from the sorted tableaux, which you'll only ever get
        via FerrersDiagram.gen_tableaux.
        '''
        self.idx = idx

    # this should be a static/class function
    @staticmethod
    def valid_static(shape, contents):
        i = 0
        rows = []
        curr_row_len = shape[0]
        prev_row_len = None
        curr_row_idx = 0
        row_idx = 0

        for i in range(len(contents)):
            if row_idx == curr_row_len:
                # start of a new row
                curr_row_idx += 1
                prev_row_len = curr_row_len
                curr_row_len = shape[curr_row_idx]
                row_idx = 0
            else:
                if row_idx < curr_row_len and i > 0:
                    if contents[i] < contents[i-1]:
                        return False
            if prev_row_len is not None:
                if contents[i] < contents[i - prev_row_len]:
                    return False
            row_idx += 1
        return True

    def valid(self):
        '''
        Checks if this young tableaux is valid.
        IE: each row must be increasing. each col must be increasing.
        '''
        for row in self.contents:
            # assert increasing
            for i in range(1, len(row)):
                if row[i-1] > row[i]:
                    return False

        max_cols = len(self.contents[0])
        for c in range(max_cols):
            for row_idx in range(1, len(self.contents)):
                row = self.contents[row_idx]
                if len(row) <= c:
                    break

                prev_row = self.contents[row_idx - 1]
                if prev_row[c] > row[c]:
                    return False
        return True

    def get_row(self, val):
        return self._row[val]

    def get_row2(self, val):
        '''
        Return the row number (1-indexed) of val
        If val is not in the tableaux (bigger than size or less than 1), return None
        '''
        for row_idx, row in enumerate(self.contents):
            if val in row:
                return row_idx + 1
        return None

    def get_col(self, val):
        return self._col[val]

    def get_col2(self, val):
        '''
        Return the row number (1-indexed) of val
        If val is not in the tableaux (bigger than size or less than 1), return None
        '''
        for c in range(len(self.contents[0])):
            col = [row[c] for row in self.contents if len(row) > c]
            if val in col:
                return c + 1
        return None

    def __repr__(self):
        '''
        Pretty prints the Young Tableau
        '''
        rep_str = ''
        for idx, l in enumerate(self.contents):
            rep_str += ''.join(map(lambda x: '[{}]'.format(x), l))
            if idx != len(self.contents) - 1:
                rep_str += '\n'
        return rep_str


    def transpose(self, transposition):
        '''
        Returns the Young Tableau you'd get by applying the transposition to this tableau if
        the resulting tableaux is valid. Otherwise, return None
        '''
        i, j = transposition
        swapped = tuple(swap(k, i, j) for k in self.vals)
        return YoungTableau.CACHE[self.partition].get(swapped, None)

    def ax_dist(self, i, j):
        '''
        Compute the Rockmore/Diaconis' axial distance for the transposition (i, j)
        in this tableaux

        i: int, element of the young tableaux
        j: int, element of the young tableaux
        '''
        row_i = self.get_row(i)
        col_i = self.get_col(i)
        row_j = self.get_row(j)
        col_j = self.get_col(j)

        return (col_j - col_i) + (row_i - row_j)

    def dist(self, i, j):
        return self.content(j) - self.content(i)

    def content(self, x):
        '''
        x: int, element of the tableaux
        This is column index minus row index of x
        '''
        col_x = self.get_col(x)
        row_x = self.get_row(x)
        return col_x - row_x

def bdown(ptup):
    bds = []
    for idx in range(len(ptup)):
        if idx < len(ptup) -1 and ptup[idx] > ptup[idx+1]:
            lst = list(ptup)
            lst[idx] = lst[idx] - 1
            bds.append(tuple(lst))
        elif idx == len(ptup) - 1:
            if ptup[idx] == 1:
                bds.append(ptup[:-1])
            else:
                lst = list(ptup)
                lst[idx] = lst[idx] - 1
                bds.append(tuple(lst))

    return bds

def test_ferrer():
    f = FerrersDiagram.from_partition((4, 2,2,1))
    print(f)
    for p in f.branch_down():
        print('---------')
        print(p)

if __name__ == '__main__':
    test_ferrer()
