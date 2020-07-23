import numpy as np
from scipy.sparse import identity, csr_matrix

from .perm import *
from .utils import hook_length, cycle_to_adj_transpositions
from .tableau import compute_syt

class SnIrrep:
    TRANS_CACHE = {}
    TRANS_SPARSE_CACHE = {}
    def __init__(self, partition, fmt='dense'):
        if partition not in SnIrrep.TRANS_CACHE:
            SnIrrep.TRANS_CACHE[partition] = {}
        if partition not in SnIrrep.TRANS_SPARSE_CACHE:
            SnIrrep.TRANS_SPARSE_CACHE[partition] = {}

        self.partition = partition
        self.n = sum(partition)
        self._tableaux = compute_syt(partition)
        self._trans_irreps = []
        self._fmt = fmt
        self._dim = int(hook_length(partition))
        self._compute_adj_trans()

        if self._fmt == 'sparse':
            self._trans_irrep = self.compute_trans_sp
        elif self._fmt == 'dense':
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
                if self._fmt == 'sparse':
                    mat = self.compute_trans_sp(*t)
                elif self._fmt == 'dense':
                    mat = self.compute_trans(*t)

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
        if self.partition in SnIrrep.TRANS_CACHE and (a, b) in SnIrrep.TRANS_CACHE[self.partition]:
            return SnIrrep.TRANS_CACHE[self.partition][(a, b)]

        size = self._dim
        mat = np.zeros((size, size))

        for tab in self.tableaux:
            other = tab.transpose(a, b)
            dist = tab.ax_dist(a, b)
            i = tab.index
            mat[i, i] = 1. / dist
            if other is not None:
                j = other.index
                if j < i:
                    continue

                mat[i, j] = np.sqrt(1 - (1. / dist) ** 2)
                mat[j, i] = mat[i, j]

        SnIrrep.TRANS_CACHE[self.partition][(a, b)] = mat
        return mat

    def compute_trans_sp(self, a, b):
        if self.partition in SnIrrep.TRANS_SPARSE_CACHE and (a, b) in SnIrrep.TRANS_SPARSE_CACHE[self.partition]:
            return SnIrrep.TRANS_SPARSE_CACHE[self.partition][(a, b)]

        rows = []
        cols = []
        data = []
        for tab in self.tableaux:
            other = tab.transpose(a, b)
            dist = tab.ax_dist(a, b)
            i = tab.index

            rows.append(i)
            cols.append(i)
            data.append(1. / dist)

            if other is not None:
                j = other.index
                if j < i:
                    # should have already computed this
                    continue

                ij_val = ji_val = np.sqrt(1 - (1. / dist) ** 2)
                #jj_val = 1. / other.ax_dist(a, b)

                rows.extend([i, j])
                cols.extend([j, i])
                data.extend([ij_val, ji_val])

        sp_mat = csr_matrix((data, (rows, cols)), shape=(self.dim, self.dim))
        SnIrrep.TRANS_SPARSE_CACHE[self.partition][(a, b)] = sp_mat
        return sp_mat

if __name__ == '__main__':
    partition = (4, 1)
    n = sum(partition)

    irrep = SnIrrep(partition, fmt='dense')
    irrep_sp = SnIrrep(partition, fmt='sparse')

    trans = [(1,2), (2,3), (3,4), (4,5)]
    for tr in trans:
        m1 = irrep.compute_trans(*tr)
        m2 = irrep_sp.compute_trans_sp(*tr).toarray()
        print(np.allclose(m1, m2))
