#!/bin/bash 

#$ -N NTD+CLR.1
#$ -pe gpu 1
#$ -q fartorgpu.q
#$ -S /bin/bash
#$ -cwd 
#$ -o NTD+CLR.1.out 
#$ -e NTD+CLR.1.err
#$ -m e
#$ -M avaldies169\@alumnes.ub.edu

# Load the modules needed 
. /etc/profile
module load gromacs/2018.1_gpu

env
export SGE_ROOT="/sge"
export CUDA_VISIBLE_DEVICES=`cat $TMPDIR/.gpus`
echo CUDA_VISIBLE_DEVICES

cd $TMPDIR
cp -r /g1819work/g19tor1/AVALDIVIA/MD/NTD+CLR/INPUTS/* $TMPDIR

#Run the job
gmx grompp -f md1.mdp -c NTD+CLR.md0.gro -p NTD+CLR.top -o NTD+CLR.md1.tpr -r NTD+CLR.md0.gro -n NTD+CLR.ndx -maxwarn 5
gmx mdrun -deffnm NTD+CLR.md1 -nt 4 

# Copy the results to our g1819work directory 
cp -r * /g1819work/g19tor1/AVALDIVIA/MD/NTD+CLR/R1/md1

wait

cd /g1819work/g19tor1/AVALDIVIA/MD/NTD+CLR/scripts/
/sge/bin/lx24-amd64/qsub md2.slm

exit
