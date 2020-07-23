import time
import pdb
import numpy as np
import random
import unittest
from perm import sn
from young_tableau import FerrersDiagram
from yor import yor
from new_yor import YorIrrep

class TestYor(unittest.TestCase):
    def test_yor_dense(self):
        partition = (4, 2, 1, 1)
        partition = (2, 1)
        s8 = sn(8)
        s8 = sn(3)

        g = random.choice(s8)
        h = random.choice(s8)
        gh = g * h

        rho = YorIrrep(partition, 'dense')
        gd = rho(g)
        hd = rho(h)
        ghd = rho(gh)
        self.assertTrue(np.allclose(gd.dot(hd), ghd))

        rho_sp = YorIrrep(partition, 'sparse')
        gsp = rho_sp(g)
        hsp = rho_sp(h)
        ghsp = rho_sp(gh)

        self.assertTrue(np.allclose(gd, gsp.toarray()))
        self.assertTrue(np.allclose(hd, hsp.toarray()))
        self.assertTrue(np.allclose(ghd, ghsp.toarray()))
        self.assertTrue(np.allclose(gsp.dot(hsp).toarray(), ghsp.toarray()))


if __name__ == '__main__':
    unittest.main()
