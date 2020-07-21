import numpy as np
import sys
import pdb
from functools import reduce
from itertools import permutations
import time

SN_CACHE = {}
SN_IDMAP = {}
HITS = {'hits': 0}
SN_TABLE = {} # n -> numpy array
SN_INV = {}
def conjugate(x, g):
    '''
    x: Perm2 object
    g: Perm2 object
    returns x conjugated by g
    '''
    return g.inv() * x * g

class ProdPerm:
    def __init__(self, *perms):
        '''
        perms: list of Perm2 objects
        '''
        self.perms = perms
        self.tup_rep = self.get_tup_rep()

    def __mul__(self, other):
        res = [a*b for a,b in zip(self.perms, other.perms)]
        return ProdPerm(*res)

    def get_tup_rep(self):
        return tuple(p.tup_rep for p in self.perms)

    def __repr__(self):
        return str(self.tup_rep)

    def inv(self):
        res = [p.inv() for p in self.perms]
        return ProdPerm(*res)

class Perm2:
    def __init__(self, p_map, n, tup_rep=None, cyc_decomp=None):
        self.size = n
        self._map = self._filled_map(p_map, self.size)
        self.cycle_decomposition = self._cycle_decomposition() if cyc_decomp is None else cyc_decomp
        self.tup_rep = self.get_tup_rep() if tup_rep is None else tup_rep
        self._id = None

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
        return Perm2(p_map, size)

    @staticmethod
    def from_tup(tup):
        if len(tup) in SN_CACHE:
            if tup in SN_CACHE[len(tup)]:
                HITS['hits'] += 1
                return SN_CACHE[len(tup)][tup]
        else:
            SN_CACHE[len(tup)] = {}

        _dict = {idx+1: val for idx, val in enumerate(tup)}
        perm = Perm2(_dict, len(tup), tup)

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
        return Perm2.from_tup(tuple(lst))

    @staticmethod
    def swap(n, a, b):
        tup = [i for i in range(1, n + 1)]
        tup[a - 1] = b
        tup[b - 1] = a
        tup = tuple(tup)
        return Perm2.from_tup(tup)

    @staticmethod
    def cont_cycle(n, a, b):
        cyc = [i for i in range(1, n+1)]
        for i in range(a, b):
            cyc[i - 1] = i + 1
        cyc[b-1] = a
        cyc = tuple(cyc)

        return Perm2.from_tup(cyc)

    @staticmethod
    def from_trans(trans, n):
        lst = [i for i in range(1, n+1)]
        lst[trans[0] - 1] = trans[1]
        lst[trans[1] - 1] = trans[0]
        tup = tuple(lst)
        return Perm2.from_tup(tup)

    @staticmethod
    def from_cycle_decomp(decomp_lst):
        n = sum(len(x) for x in decomp_lst)
        output_lst = [i for i in range(1, n+1)]
        for cyc in decomp_lst:
            for idx, val in enumerate(cyc):
                nxt_val = cyc[(idx + 1) % len(cyc)]
                output_lst[val - 1] = nxt_val
        return Perm2.from_tup(tuple(output_lst))

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
        return Perm2.from_tup(new_tup)

    # this is actually slower than the other __mul__
    #def __mul__2(self, other):
    #    prod_id = SN_TABLE[self.size][self.id, other.id] # this maps to an id
    #    return SN_IDMAP[self.size][prod_id]

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

    #def inv2(self):
    #    perm_id = SN_INV[self.size][self.id]
    #    return SN_IDMAP[self.size][perm_id]

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
            return Perm2(rev_map, self.size)
        rev_map = {v: k for k, v in self._map.items()}
        return Perm2(rev_map, self.size)

    @property
    def id(self):
        return self._id

    def set_id(self, _id):
        self._id = _id
 
def sn(n, prefix='/local/hopan/'):
    # load mult table?
    if n in SN_CACHE and len(SN_CACHE[n]) == np.math.factorial(n):
        return list(SN_CACHE[n].values())

    perm_tups = permutations(range(1, n+1))
    perms = [Perm2.from_tup(t) for t in perm_tups]
    #SN_IDMAP[n] = {}
    for idx, p in enumerate(perms):
        p.set_id(idx)
        #SN_IDMAP[n][idx] = p

    SN_CACHE[n] = {p.tup_rep: p for p in perms}
    #print('using mul cache')
    #print('loading: {}'.format(prefix + 's{}_table.npy'.format(n)))
    #SN_TABLE[n] = np.load(prefix + 's{}_table.npy'.format(n))
    #SN_INV[n] = np.load(prefix + 's{}_inv.npy'.format(n))
    #print('done loading: {}'.format(prefix + 's{}_table.npy'.format(n)))

    return perms

'''
def test():
    for n in range(2,11):
        print('n = {}'.format(n))
        f = lambda x, y: x*y
        start = time.time()
        for i in range(10000):
            ps = [Perm([(i, i+1)]) for i in range(1, n+1)]
            p = reduce(f, ps)
        end = time.time()
        print('Time for orginal perm: {:.2f}'.format(end - start))

        start =time.time()
        for i in range(10000):
            ps = [Perm2({i:i+1, i+1:i}, n) for i in range(1, n+1)]
            p = reduce(f, ps)
        end = time.time()
        print('Time for perm with maps: {:.2f}'.format(end - start))

        print('=' * 10)
    #start = time.time()
    #for i in range(10000):
    #    ps = [Perm([(i, i+1)]) for i in range(1, n+1)]
    #    p = reduce(f, ps)
    #end = time.time()
    #print('Time for orginal perm 2nd time: {:.2f}'.format(end - start))
'''

def mult_table(n, inv_save, table_save):
    _sn = sn(n)
    table = np.zeros((len(_sn), len(_sn)), dtype=np.uint16)
    inv_table = np.zeros(len(_sn), dtype=np.uint16)
    for idx, p in enumerate(_sn):
        p.set_id(idx)

    for i, p in enumerate(_sn):
        inv_table[i] = p.inv().id

        for j, k in enumerate(_sn):
            table[i, j] = (p * k).id
            table[j, i] = (k * p).id

    np.save(inv_save.format(n), inv_table)
    np.save(table_save, table)

if __name__ == '__main__':
    n = int(sys.argv[1])
    start = time.time()
    mult_table(n)
    end = time.time()
    print('s{} | Elapsed: {:.2f}'.format(n, end - start))
