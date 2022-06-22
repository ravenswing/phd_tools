import os
import subprocess
import time
import sys
import glob

#svr = 'zcqsrev@thomas.rc.ucl.ac.uk'
svr = 'rxe15-ffg01@jade.hartree.stfc.ac.uk'
svrdir = 'baker_peptides/CHARMM36m/'
srcdir = 'baker_peptides/CHARMM36m/'
mols = ['KKKK', 'KKEE']

watermodel = [
        #"amber14sb_tip4pd.ff/"
        #"amber_disp.ff/"
        "charmm36m_sw.ff/"
       ]

ext = False

#-----------------------------------------------------------------------------

for mol in mols:
    # make directory
    try:
        subprocess.call(
            'ssh {} mkdir -p {}:~/{}{}'.format(svr,svrdir,mol),
        shell=True)
    except:
        print( 'Unable to make directory: {}'.format(mol))

    # copy files to new directory
    if ext:
        try: 
            subprocess.call(
            'scp {}{m}/mdrun.tpr {}:~/{}{m}/'.format(srcdir,svr,svdir,m=mol),
            shell=True)
            
            paths = glob.glob('./{}{}/*'.format(srcdir,mol))
            print( paths )
            select = [ line[:-4].split('/') for line in paths if line[0] != '#' and '.edr' in line]
            print( select )
            stem =  sorted(select,reverse=True)[0][-1]
            print (stem)

            subprocess.call(
            'scp {}{m}/{}.cpt {}:~/{}{m}/'.format(srcdir,stem,svr,m=mol),
            shell=True)
        except:
            print( 'Unable to transfer top/gro files')

        # edit working directory in .sh run file
        with open('Scripts/MN/mdrun_ex.sh', 'r') as file: 
            data = file.readlines()
        #  print data
        print( 'Editing line:  {}'.format(data[6]))
        nxt = stem[:-1]+str( int( stem[-1] ) +1) 
        data[1] = '#SBATCH --job-name="MD{}_{}" \n'.format(nxt[-1], mol)
        data[13] = 'mpirun $GMX  mdrun -v -s MDrun_ex.tpr -cpi {}.cpt -deffnm {} -append -maxh 23.5 \n'.format(stem, nxt) 
        with open('Scripts/MN/mdrun_ex.sh', 'w') as file:
            file.writelines( data )

    else:
        try:
            subprocess.call(
            'scp {}{m}/05-MD/{m}_1barframe.gro {}:~/{}{m}/'.format(srcdir,svr,svrdir,m=mol),
            shell=True)
            subprocess.call(
            'scp {}{m}/05-MD/{m}.top {}:~/{}{m}/'.format(srcdir,svr,svrdir,m=mol),
            shell=True)
            subprocess.call(
            'scp -r {}{m}/05-MD/{} {}:~/{}{m}/'.format(srcdir,watermodel[0],svr,svrdir,m=mol),
            shell=True)
            subprocess.call(
            'scp {}{m}/05-MD/min.tpr {}:~/{}{m}/'.format(srcdir,svr,svrdir,m=mol),
            shell=True)
            subprocess.call(
            'scp {}{m}/00-Prep/{m}_initial.pdb {}:~/{}{m}/'.format(srcdir,svr,svrdir,m=mol),
            shell=True)
            subprocess.call(
            'scp -r {}{m}/05B-SWISH/* {}:~/{}{m}/'.format(srcdir,svr,svrdir,m=mol),
            shell=True)
        except:
            print( 'Unable to transfer top/gro files')
            # edit working directory in .sh run file

        # edit working directory in .sh run file
        for fl in ['mdrun.sh', 'mdcont.sh']:
            with open('Scripts/TH/{}'.format(fl), 'r') as file: 
                data = file.readlines()
            #  print data
            print( 'Editing line:  {}'.format(data[6])) 
            data[5] = '#$ -N MD_{} \n'.format(mol)
            data[8] = '#$ -wd /home/zcqsrev/{}{} \n'.format(svrdir,mol)
            with open('Scripts/TH/{}'.format(fl), 'w') as file:
                file.writelines( data )
    # copy (EDITED) .sh and .mdp files to Grace
    try:
        subprocess.call(
        'scp Scripts/TH/* {}:~/{}{}/'.format(svr,svrdir,mol),
        shell=True)
    except:
        print( 'Unable to transfer mdrun files')

    print( '{}  successfully prepared \n'.format(mol))

print( '----------  All {} molecules successfully prepared  ----------'.format(len(mols)))

print( '\nStarting batch submission... \n')
for mol in mols:

    sub = subprocess.Popen( 'ssh -T {}'.format(svr),
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True, bufsize=0,
        shell=True )
  
    sub.stdin.write('ls .\n')
    sub.stdin.write('cd {}{}/\n'.format(svrdir,mol))
    sub.stdin.write('pwd\n')
#    sub.stdin.write('qsub {fn}.sh\n'.format(fn='mdrun_ex' if ext else 'mdrun'))
    out, errors = sub.communicate()
    print (out)
    print (errors)
              
