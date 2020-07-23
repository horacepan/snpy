# snpy
SnPy is a simple library for computing irreducible representation matrices (using Young's orthogonal representation) of the the symmetric group.

## Requirements
- numpy
- scipy

Sample Usage
```
from snpy.perm import Perm
from snpy.sn_irrep import SnIrrep

partition = (4, 1)
rho = SnIrrep(partition, fmt='dense')
g = Perm.cont_cycle(1, 3, 5)
g_rep = rho(g)
```
