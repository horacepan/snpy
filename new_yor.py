import numpy as np
from YOR import yor
from young_tableau import FerrersDiagram
from scipy.sparse import identity
from perm import *
from utils import hook_length, cycle_to_adj_transpositions
from tableau import compute_syt

YOR_TRANS_CACHE = {}
class YOR:
    TRANS_CACHE = {}
    def __init__(self, partition, fmt='dense'):
        if partition not in YOR_TRANS_CACHE:
            YOR_TRANS_CACHE[partition] = {}

        self.partition = partition
        self.n = sum(partition)
        self._tableaux = compute_syt(partition)
        self._trans_irreps = []
        self._fmt = fmt
        self._dim = int(hook_length(partition))
        self._compute_adj_trans()

        if self._fmt == 'sparse':
            self._trans_irrep = self.compute_trans_sp
        elif self._fmt == 'sparse':
            self._trans_irrep = self.compute_trans
        else:
            pass

    @property
    def dim(self):
        return self._dim

    @property
    def tableaux(self):
        return self._tableaux

    def __call__(self, perm):
        # do the bubble sort multiplication
        if self._fmt == 'dense':
            irrep = np.eye(self._dim)
        elif self._fmt == 'sparse':
            irrep = identity(self._dim, format='csr')

        for cycle in perm.cycle_decomposition:
            for t in cycle_to_adj_transpositions(cycle, self.n):
                mat = self.compute_trans_sp(*t) if self._fmt == 'sparse' else self.compute_trans(*t)
                irrep = irrep.dot(mat)

        return irrep

    def _compute_adj_trans(self):
        adj_trans = [(i, i+1) for i in range(1, self.n)]

        for (a, b) in adj_trans:
            if self._fmt == 'dense':
                mat = self.compute_trans(a, b)
            elif self._fmt == 'sparse':
                mat = self.compute_trans_sp(a, b)
            else:
                raise Exception('Only allowed formats are sparse(csr) and dense(np).')

    def compute_trans(self, a, b):
        if self.partition in YOR_TRANS_CACHE and (a, b) in YOR_TRANS_CACHE[partition]:
            return YOR_TRANS_CACHE[self.partition][(a, b)]

        size = self._dim
        rep = np.zeros((size, size))

        for tab in self.tableaux:
            other = tab.transpose(a, b)
            dist = tab.ax_dist(a, b)
            i = tab.index
            rep[i, i] = 1. / dist

            if other is not None:
                j = other.index
                rep[i, j] = np.sqrt(1 - (1. / dist) ** 2)
                rep[j, i] = rep[i, j]
                rep[j, j] = 1. / other.ax_dist(a, b)

        try:
            YOR_TRANS_CACHE[self.partition][(a, b)] = rep
        except:
            pdb.set_trace()
        return rep

    def compute_trans_sp(self, a, b):
        if self.partition in YOR_TRANS_CACHE and (a, b) in YOR_TRANS_CACHE[partition]:
            return YOR_TRANS_CACHE[self.partition][(a, b)]

        rows = []
        cols = []
        data = []
        for tab in self.tableaux:
            other = tab.transpose(transposition)
            dist = tab.ax_dist(*transposition)
            i = tab.index

            rows.append(i)
            cols.append(i)
            data.append(1. / dist)

            if other is not None:
                j = other.index

                ij_val = ji_val = np.sqrt(1 - (1. / dist) ** 2)
                jj_val = 1. / other.ax_dist(*transposition)

                rows.extend([i, j, j])
                cols.extend([j, i, j])
                data.extend([ij_val, ji_val, jj_val])

        sp_mat = csr_matrix(data, (rows, cols), shape=(self.dim, self.dim))
        return sp_mat

if __name__ == '__main__':
    partition = (4, 1)
    n = sum(partition)

    irrep = YOR(partition, fmt='dense')
    perm = Perm.cycle(1, 4, n)
    m1 = irrep(perm)
    m2 = irrep(perm.inv())
    y1 = yor(FerrersDiagram(partition), perm)
    y2 = yor(FerrersDiagram(partition), perm.inv())
    print(np.allclose(m1, y1))
    print(np.allclose(m2, y2))
