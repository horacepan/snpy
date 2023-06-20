## snpy
SnPy is a simple library for computing irreducible representation matrices (using Young's orthogonal representation) of the symmetric group.

## Requirements
- numpy
- scipy

## Installation
```
python setup.py install
```

## Sample Usage
```
from snpy.perm import Perm
from snpy.sn_irrep import SnIrrep

partition = (4, 1)
rho = SnIrrep(partition, fmt='dense')
g = Perm.cycle(1, 3, 5)
print(rho(g))
```
See the sample [Jupyter notebook](https://github.com/horacepan/snpy/blob/master/example.ipynb) for more details on how to generate various permutations and their corresponding irrep matrices.
