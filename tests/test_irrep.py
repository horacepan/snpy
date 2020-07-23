import time
import random
import unittest

import numpy as np
from snpy.perm import Perm, sn
from snpy.sn_irrep import SnIrrep

class TestSnIrrep(unittest.TestCase):
    def test_yor_dense(self):
        partition = (4, 2, 1, 1)
        s8 = sn(8)

        g = random.choice(s8)
        h = random.choice(s8)
        gh = g * h

        rho = SnIrrep(partition, 'dense')
        gd = rho(g)
        hd = rho(h)
        ghd = rho(gh)
        self.assertTrue(np.allclose(gd.dot(hd), ghd))

        rho_sp = SnIrrep(partition, 'sparse')
        gsp = rho_sp(g)
        hsp = rho_sp(h)
        ghsp = rho_sp(gh)

        self.assertTrue(np.allclose(gd, gsp.toarray()))
        self.assertTrue(np.allclose(hd, hsp.toarray()))
        self.assertTrue(np.allclose(ghd, ghsp.toarray()))
        self.assertTrue(np.allclose(gsp.dot(hsp).toarray(), ghsp.toarray()))

    def test_yor_inv(self):
        partition = (4, 3, 1)
        s8 = sn(8)

        rho = SnIrrep(partition, 'dense')
        g = random.choice(s8)
        gmat = rho(g)
        mat_inv = np.linalg.inv(gmat)
        ginv_mat = rho(g.inv())

        self.assertTrue(np.allclose(mat_inv, ginv_mat))

    def test_s3(self):
        partition = (3, 1)
        rho = SnIrrep(partition, fmt='dense')
        p12 = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
        p23 = np.array([[0.5, np.sqrt(3)/2., 0], [np.sqrt(3)/2., -0.5, 0], [0, 0, 1]])
        p34 = np.array([[1, 0, 0], [0, 1./3., np.sqrt(8)/3], [0, np.sqrt(8)/3, -1./3.]])

        self.assertTrue(np.allclose(p12, rho(Perm.trans(1, 2, 4))))
        self.assertTrue(np.allclose(p23, rho(Perm.trans(2, 3, 4))))
        self.assertTrue(np.allclose(p34, rho(Perm.trans(3, 4, 4))))

if __name__ == '__main__':
    unittest.main()
