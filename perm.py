import numpy as np
import sys
import pdb
from functools import reduce
from itertools import permutations
import time

class Perm:
    def __init__(self, lst_map):
        self.size = len(lst_map)
        self._lst_map = lst_map
        self._tup_rep = tuple(self._lst_map)
        self._cycle_decomposition = self.compute_cycle_decomposition()

    @staticmethod
    def eye(size):
        return Perm([i for i in range(1, size + 1)])

    @staticmethod
    def from_cycle_decomp(decomp_lst, n):
        lst = [i for i in range(1, n+1)]
        for cyc in decomp_lst:
            for idx, val in enumerate(cyc):
                nxt_val = cyc[(idx + 1) % len(cyc)]
                lst[val - 1] = nxt_val
        return Perm(lst)

    @staticmethod
    def cont_cycle(st, end, n):
        cyc = [i for i in range(1, n+1)]
        for i in range(st, end):
            cyc[i - 1] = i + 1
        cyc[b-1] = a

        return Perm(cyc)

    @staticmethod
    def trans(a, b, n):
        lst = [i for i in range(1, n+1)]
        lst[a - 1] = b
        lst[b - 1] = a
        return Perm(lst)

    @property
    def tup(self):
        return self._tup_rep

    @property
    def cycle_decomposition(self):
        return self._cycle_decomposition

    def inv(self):
        inv_map = [self._lst_map.index(i + 1) + 1 for i in range(self.size)]
        return Perm(inv_map)

    def mat(self):
        pmat = np.zeros((self.size, self.size))
        for j in range(self.size):
            i = self._lst_map[j] - 1
            pmat[i, j] = 1

        return pmat

    def mul(self, perm):
        return self * perm

    def __call__(self, x):
        return self.tup[x-1]

    def __getitem__(self, x):
        return self.tup[x-1]

    def __repr__(self):
        return str(self.cycle_decomposition)

    def __mul__(self, other):
        if self.size != other.size:
            raise Exception('Currently cant mult two perms of diff sizes!')

        g = self._lst_map
        h = other._lst_map
        new_lst = [g[h[i] - 1] for i in range(self.size)]
        return Perm(new_lst)

    def __len__(self):
        return self.size

    def __eq__(self, other):
        return isinstance(self, type(other)) and self._lst_map == other._lst_map

    def compute_cycle_decomposition(self):
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
                    curr = self._lst_map[curr-1]
                else:
                    if len(curr_cycle) > 1:
                        cyc_decomp.append(curr_cycle)
                    curr_cycle = []
                    break
            seen.add(curr)

        return cyc_decomp

def sn(n):
    perm_tups = permutations(range(1, n+1))
    return [Perm(p) for p in perm_tups]

if __name__ == '__main__':
    x = Perm.from_cycle_decomp([(1,2), (3, 4, 5)])
    y = Perm.trans(1, 2, 5)
    print(x)
    print(x.inv())
    print(y)
    print(x * y)
