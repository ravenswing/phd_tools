import os
import subprocess
import time
import sys

mn = 'pr1eeq09@mn3.bsc.es' 

# set molecule name from user input
if len(sys.argv) < 2:
    print("\nPlease try again with valid molecule name...\n")
    sys.exit()
mol = sys.argv[1]
print ('Attempting to restart:  {}'.format(mol))
# sync working folder
try:
    subprocess.call(
        'rsync -avzhPe ssh --exclude "*.trr" {}:/gpfs/scratch/pr1eeq00/pr1eeq09/{m}/ Hydrophobic/{m}/06-RunMD/'.format(mn,m=mol),
        shell=True)
except:
    print( 'Unable to make sync MareNostrum directory for {}'.format(mol))
# select the next run number
ls = subprocess.Popen(['ssh','{}'.format(mn), 'ls /gpfs/scratch/pr1eeq00/pr1eeq09/{}/'.format(mol)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err =  ls.communicate()
output = out.decode('utf-8').splitlines()
select = [ line[:-4] for line in output if line[0] != '#' and '.edr' in line]
if len(select) < 1 or len(select) == 1 and select[0] != 'MDrun':
    sys.exit('Unable to restart: First MD run output not found')
if len(select)==1 and select[0]=='MDrun': 
    newline = 'MDrun01'
else:
    newline = 'MDrun0{}'.format(int( sorted(select,reverse=True)[0][-1] )+1)
print ('Next run:  {}'.format(newline))
# edit the restart script to have correct working directory and run number
nc = 'gerun $GMX  mdrun -v -s mdrun.tpr -deffnm {} -tunepme -pin on -cpi {}.cpt  -append -maxh 23.5\n'.format(newline, sorted(select,reverse=True)[0])
# edit working directory in .sh run file
with open('Scripts/Grace/mdrst.sh', 'r') as file:
    data = file.readlines()
print(data)
print( 'Editing line:  {}'.format(data[24]))
data[1] ='#$ -N MD{}_{} \n'.format(newline[-1],mol)
data[6] = '#$ -wd /home/zcqsrev/Scratch/{} \n'.format(mol)
data[24] = nc 
with open('Scripts/Grace/mdrst.sh', 'w') as file:
    file.writelines( data )
# copy the EDITED rst script to Grace
try:
    subprocess.call(
        'scp Scripts/Grace/mdrst.sh {}:~//gpfs/scratch/pr1eeq00/pr1eeq09/{}/'.format(mn,mol),
        shell=True)
except:
    print ('Unable to transfer mdrun files')
print ('\nStarting batch submission... \n')
# submit the new job
sub = subprocess.Popen( 'ssh -T {}'.format(mn),
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    universal_newlines=True, bufsize=0,
    shell=True )
sub.stdin.write('cd /gpfs/scratch/pr1eeq00/pr1eeq09/{}/\n'.format(mol))
sub.stdin.write('pwd\n')
sub.stdin.write('qsub mdrst.sh\n')
out, errors = sub.communicate()
print (out)
print (errors)

