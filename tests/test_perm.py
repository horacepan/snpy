import random
import unittest

#from .context import snpy
from snpy.perm import Perm, sn

class TestPerm(unittest.TestCase):
    def test_eye(self):
        n = 10
        g = Perm.eye(n)
        self.assertEqual(g.tup, tuple(range(1, n + 1)))
        for i in range(1, n + 1):
            self.assertEqual(g[i], i)

    def test_trans(self):
        g = Perm.trans(1, 4, 6)
        self.assertEqual(g.tup, (4, 2, 3, 1, 5, 6))

    def test_cycle(self):
        g = Perm.cycle(1, 5, 5)
        self.assertEqual(g.tup, (2, 3, 4, 5, 1))

    def test_mul(self):
        g = Perm.trans(1, 2, 6)
        h = Perm.trans(3, 4, 6)
        gh = g * h
        self.assertEqual(gh.tup, (2, 1, 4, 3, 5, 6))

    def test_inv(self):
        g = Perm((3, 2, 4, 1, 5))
        ginv = Perm((4, 2, 1, 3, 5))

        self.assertEqual(g.inv().tup, ginv.tup)

    def test_cycle_decomp(self):
        g = Perm.cycle(1, 3, 10)
        h = Perm.trans(5, 7, 10)
        f = Perm.trans(4, 6, 10)

        ghf = g * h * f
        self.assertEqual(ghf.cycle_decomposition, [[1, 2, 3], [4, 6], [5, 7]])

    def test_mat(self):
        n = 8
        perms = sn(n)
        g = random.choice(perms)
        mat = g.mat()

        for j in range(1, n + 1):
            i = g[j]
            self.assertEqual(mat[i-1, j-1], 1)

if __name__ == '__main__':
    unittest.main()
