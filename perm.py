import numpy as np
import sys
import pdb
from functools import reduce
from itertools import permutations
import time

SN_CACHE = {}
HITS = {'hits': 0}
SN_TABLE = {} # n -> numpy array
SN_INV = {}

class Perm:
    def __init__(self, p_map, n, tup_rep=None, cyc_decomp=None):
        self.size = n
        self._map = self._filled_map(p_map, self.size)
        self.cycle_decomposition = self._cycle_decomposition() if cyc_decomp is None else cyc_decomp
        self.tup_rep = self.get_tup_rep() if tup_rep is None else tup_rep

        # add permutation to the cache
        if self.size in SN_CACHE:
            if (self.tup_rep not in SN_CACHE):
                SN_CACHE[self.size][self.tup_rep] = self
        elif self.size not in SN_CACHE:
            SN_CACHE[self.size] = {}
            SN_CACHE[self.size][self.tup_rep] = self

    def _filled_map(self, _map, size):
        for i in range(1, size+1):
            if i not in _map:
                _map[i] = i
        return _map

    def get_tup_rep(self):
        lst = []
        for i in range(1, len(self._map) + 1):
            lst.append(self._map.get(i, i))
        return tuple(lst)

    @staticmethod
    def eye(size):
        p_map = {i:i for i in range(1, size+1)}
        return Perm(p_map, size)

    @staticmethod
    def from_tup(tup):
        if len(tup) in SN_CACHE:
            if tup in SN_CACHE[len(tup)]:
                HITS['hits'] += 1
                return SN_CACHE[len(tup)][tup]
        else:
            SN_CACHE[len(tup)] = {}

        _dict = {idx+1: val for idx, val in enumerate(tup)}
        perm = Perm(_dict, len(tup), tup)

        # store the thing
        SN_CACHE[len(tup)][perm.tup_rep] = perm
        return perm

    @staticmethod
    def from_cycle_decomp(decomp_lst):
        n = sum(len(x) for x in decomp_lst)
        lst = [i for i in range(1, n+1)]
        for cyc in decomp_lst:
            for idx, val in enumerate(cyc):
                nxt_val = cyc[(idx + 1) % len(cyc)]
                lst[val - 1] = nxt_val
        return Perm.from_tup(tuple(lst))

    @staticmethod
    def swap(n, a, b):
        tup = [i for i in range(1, n + 1)]
        tup[a - 1] = b
        tup[b - 1] = a
        tup = tuple(tup)
        return Perm.from_tup(tup)

    @staticmethod
    def cont_cycle(n, a, b):
        cyc = [i for i in range(1, n+1)]
        for i in range(a, b):
            cyc[i - 1] = i + 1
        cyc[b-1] = a
        cyc = tuple(cyc)

        return Perm.from_tup(cyc)

    @staticmethod
    def from_trans(trans, n):
        lst = [i for i in range(1, n+1)]
        lst[trans[0] - 1] = trans[1]
        lst[trans[1] - 1] = trans[0]
        tup = tuple(lst)
        return Perm.from_tup(tup)

    @staticmethod
    def from_cycle_decomp(decomp_lst):
        n = sum(len(x) for x in decomp_lst)
        output_lst = [i for i in range(1, n+1)]
        for cyc in decomp_lst:
            for idx, val in enumerate(cyc):
                nxt_val = cyc[(idx + 1) % len(cyc)]
                output_lst[val - 1] = nxt_val
        return Perm.from_tup(tuple(output_lst))

    def __call__(self, x):
        return self._map.get(x, x)

    def __getitem__(self, x):
        return self._map.get(x, x)

    def __repr__(self):
        # mostly for debugging so dont care about recomputing this conversion
        return str(self.cycle_decomposition)

    def __mul__(self, other):
        g = self.tup_rep
        h = other.tup_rep
        if self.size != other.size:
            raise Exception('Currently cant mult two perms of diff sizes!')

        new_tup = tuple(g[h[i] - 1] for i in range(self.size))
        return Perm.from_tup(new_tup)

    def __len__(self):
        return self.size

    def __hash__(self):
        return hash(self.tup_rep)

    def __eq__(self, other):
        return isinstance(self, type(other)) and self.tup_rep == other.tup_rep

    def _cycle_decomposition(self):
        cyc_decomp = []
        curr_cycle = []
        seen = set()
        for i in range(1, self.size+1):
            if i in seen:
                continue
            curr = i

            while True:
                if curr not in seen:
                    curr_cycle.append(curr)
                    seen.add(curr)
                    curr = self._map.get(curr, curr)
                else:
                    if len(curr_cycle) > 1:
                        cyc_decomp.append(curr_cycle)
                    curr_cycle = []
                    break
            seen.add(curr)

        return cyc_decomp

    def to_tup(self):
        return tuple(self._map[i] for i in range(1, self.size+1))

    def inv(self):
        rev_lst = [0] * self.size
        for idx, v in enumerate(self.tup_rep):
            # idx is 0 indexed, but v is 1 indexed?
            rev_lst[v-1] = idx + 1
        rev_tup = tuple(rev_lst) # surely theres a better way to do this?

        if rev_tup in SN_CACHE[self.size]:
            HITS['hits'] += 1
            return SN_CACHE[self.size][rev_tup]
        else:
            rev_map = {v: k for k, v in self._map.items()}
            return Perm(rev_map, self.size)
        rev_map = {v: k for k, v in self._map.items()}
        return Perm(rev_map, self.size)

    def mat(self):
        pmat = np.zeros((self.size, self.size))
        for j in range(self.size):
            i = self.tup_rep[j] - 1
            pmat[i, j] = 1

        return pmat

def sn(n):
    perm_tups = permutations(range(1, n+1))
    return [Perm.from_tup(p) for p in perm_tups]

if __name__ == '__main__':
    s4 = sn(4)
    print(s4)
    x = Perm.from_tup((2,3,4,5,1))
    print(x.mat())
