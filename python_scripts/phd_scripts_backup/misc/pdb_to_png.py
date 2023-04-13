import sys
import __main__
import pymol


srcdir = 'Hydrophilic/KKEE/06-RunMD/'
cutoff = '0.6'
method = 'gr'
'''
#for i in range(2,6):
def createPNG(cluster_num):

    nm = 'clst_{}_{}.pdb0{}'.format(cutoff, method, cluster_num)
    print(nm)


    pdb_file = srcdir+nm+'.pdb'
    pdb_name = nm
    print('.')
    print(pdb_file, pdb_name)
    pymol.cmd.load(pdb_file, pdb_name)
    pymol.cmd.disable("all")
    print('.')
    pymol.cmd.enable(pdb_name)
    print( pymol.cmd.get_names())
    pymol.cmd.hide('all')
    pymol.cmd.show('cartoon', str(pdb_name))
    pymol.cmd.set('ray_opaque_background', 0)
    pymol.cmd.color('red', 'ss h')
    pymol.cmd.color('yellow', 'ss s')
    pymol.cmd.png("{}{}.png".format(srcdir,pdb_name))

sys.path.append('/c5/shared/pymol/1.7.0.0-python-2.7.5-shared/lib/python2.7/site-packages/')
__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
pymol.finish_launching()
print('pymol lauched')
for i in range(1,6):
    createPNG(i)

pymol.cmd.quit()
'''

def createPNG(cluster_num):

    sys.path.append('/c5/shared/pymol/1.7.0.0-python-2.7.5-shared/lib/python2.7/site-packages/')
    nm = 'clst_{}_{}.pdb0{}'.format(cutoff, method, cluster_num)
    print(nm)

    __main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
    pymol.finish_launching()
    print('pymol lauched')

    pdb_file = srcdir+nm+'.pdb'
    pdb_name = nm
    print('.')
    print(pdb_file, pdb_name)
    pymol.cmd.load(pdb_file, pdb_name)
    pymol.cmd.disable("all")
    print('.')
    pymol.cmd.enable(pdb_name)
    print( pymol.cmd.get_names())
    pymol.cmd.hide('all')
    pymol.cmd.show('cartoon')
    pymol.cmd.set('ray_opaque_background', 0)
    pymol.cmd.color('red', 'ss h')
    pymol.cmd.color('yellow', 'ss s')
    pymol.cmd.png("{}{}.png".format(srcdir,pdb_name))
    pymol.cmd.quit()


createPNG(3)
