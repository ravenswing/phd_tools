import subprocess
import glob
from itertools import chain

gro_name = 'FRG_GMX.gro'

path = glob.glob('./*built*')
fname =  [s for s in path if '#' not in s][0].split('/')[-1]
mol = fname.split('_')[0]


with open('./{}'.format(gro_name),'r') as f:
    ligand_gro = f.readlines()
with open('./{}'.format(fname),'r') as f:
    protein_gro = f.readlines()

total = (int(ligand_gro[1].split('\n')[0]) + int(protein_gro[1].split('\n')[0]))
line2 = ' '+str(total)+'\n'
both_gro = chain( protein_gro[0], line2, protein_gro[2:-1],ligand_gro[2:-1], protein_gro[-1])
with open('./{}_edited.gro'.format(mol), 'w+') as f:
   f.writelines(str(l) for l in both_gro)

with open('./topol.top') as f:
    d = '['
    topol = [d+e for e in f.read().split(d) if e]
topol[0] = topol[0][1:]
include1 = '; Include ligand atom types \n#include \"./ligand_atomtypes.itp\" \n\n'
include2 = '\n; Include ligand topology   \n#include \"./ligand.itp\" \n\n'
include3 = 'FRG              1\n'
n1 = [i for i,s in enumerate(topol) if 'moleculetype' in s][0]
n2 = [i for i,s in enumerate(topol) if 'posre.itp' in s][0]
split = topol[n2].split('#endif')
print('n2')
print(topol[n2])
print('n2+1')
print(topol[n2+1])
print('n2+2')
print(topol[n2+2])
new_top = chain(topol[0:n1], include1, topol[n1:n2], split[0],'#endif', include2, split[1], topol[n2+1:], include3)
with open('./topol.top', 'w') as f:
    f.writelines(str(l) for l in new_top)
