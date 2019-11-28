import pandas as pd
import numpy as np
from itertools import repeat

npep = 3
nres = 4
tmp = ['pep'+str(x+1) for x in range(npep)]
pep = [x for item in tmp for x in repeat(item, nres)]

tmp = ['res'+str(x+1) for x in range(nres)]
res = []
for i in np.arange(nres):
    res.extend(tmp)
    i += 1

print(res)

arrays = [pep, res]

tuples = list(zip(*arrays))
print(tuples)

index = pd.MultiIndex.from_tuples(tuples, names=['Peptide', 'Residue'])

df = pd.DataFrame(np.random.randn(3, 12), columns=index)

print(type(df['pep1', 'res1']))

