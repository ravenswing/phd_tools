import os
import subprocess
import time
import sys
import glob

svr = 'pr1eeq09@mn3.bsc.es'
srcdir = 'Results/Hydrophobic/'
mols = ['HHHP', 'HHHP+']

ext = True

for mol in mols:
    # make directory
    try:
        subprocess.call(
        'ssh {} mkdir -p /gpfs/scratch/pr1eeq00/pr1eeq09/{}'.format(svr,mol),
        shell=True)
    except:
        print( 'Unable to make directory: {}'.format(mol))

    # copy files to new directory
    if ext:
        try: 
            subprocess.call(
            'scp {}{m}/mdrun.tpr {}:/gpfs/scratch/pr1eeq00/pr1eeq09/{m}/'.format(srcdir,svr,m=mol),
            shell=True)
            
            paths = glob.glob('./{}{}/*'.format(srcdir,mol))
            print( paths )
            select = [ line[:-4].split('/') for line in paths if line[0] != '#' and '.edr' in line]
            print( select )
            stem =  sorted(select,reverse=True)[0][-1]
            print (stem)

            subprocess.call(
            'scp {}{m}/{}.cpt {}:/gpfs/scratch/pr1eeq00/pr1eeq09/{m}/'.format(srcdir,stem,svr,m=mol),
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
            'scp {}{m}/05-MD/{m}_1barframe.gro {}:/gpfs/scratch/pr1eeq00/pr1eeq09/{m}/'.format(srcdir,svr,m=mol),
            shell=True)
            subprocess.call(
            'scp {}{m}/05-MD/{m}.top {}:/gpfs/scratch/pr1eeq00/pr1eeq09/{m}/'.format(srcdir,svr,m=mol),
            shell=True)
        except:
            print( 'Unable to transfer top/gro files')
            # edit working directory in .sh run file

        # edit working directory in .sh run file
        with open('Scripts/MN/mdrun.sh', 'r') as file: 
            data = file.readlines()
        #  print data
        print( 'Editing line:  {}'.format(data[6]))
        nxt = stem[:-1]+str( int( stem[-1] ) +1) 
        data[1] = '#SBATCH --job-name="MD{}_{}" \n'.format(nxt[-1], mol)
        data[13] = 'mpirun $GMX  mdrun -v -s MDrun_ex.tpr -cpi {}.cpt -deffnm {} -append -maxh 23.5 \n'.format(stem, nxt) 
        with open('Scripts/MN/mdrun.sh', 'w') as file:
            file.writelines( data )

        
  # edit working directory in .sh restart script
#  with open('Scripts/Grace/mdrst.sh', 'r') as file: 
 #   data = file.readlines()

  #print ('Editing line:  {}'.format(data[6]))
 # data[6] = '#$ -wd /home/zcqsrev//gpfs/scratch/pr1eeq00/pr1eeq09/{} \n'.format(mol) 
 # with open('Scripts/Grace/mdrst.sh', 'w') as file:
  #  file.writelines( data )
  
  # copy (EDITED) .sh and .mdp files to Grace
    try:
        subprocess.call(
        'scp Scripts/MN/* {}:/gpfs/scratch/pr1eeq00/pr1eeq09/{}/'.format(svr,mol),
        shell=True)
  except:
        print( 'Unable to transfer mdrun files')

  print( '{}  successfully prepared \n'.format(mol))

print( '----------  All {} molecules successfully prepared  ----------'.format(len(mols)))
# submit the first mdrun to the queue

print( '\nStarting batch submission... \n')
for mol in mols:

    sub = subprocess.Popen( 'ssh -T {}'.format(svr),
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True, bufsize=0,
        shell=True )
  
    sub.stdin.write('ls .\n')
    sub.stdin.write('cd /gpfs/scratch/pr1eeq00/pr1eeq09/{}/\n'.format(mol))
    sub.stdin.write('pwd\n')
    sub.stdin.write('sbatch {fn}.sh\n'.format(fn='mdrun_ex' if ext else 'mdrun'))
    out, errors = sub.communicate()
    print (out)
    print (errors)
              
