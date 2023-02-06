#!/bin/csh 

#$ -S /bin/bash
#$ -N BASE_NAME
#$ -pe smp 1
#$ -q iqtc10_g1819_gpu.q
#$ -cwd 
#$ -o BASE_NAME_equil.out
#$ -e BASE_NAME_equil.err

# Load the modules needed 
source /etc/profile

module load amber/20_multigpu

export OMP_NUM_THREADS=1

# Copy inputs and files needed to the directory where the jobs will run
echo " SGE_O_WORKDIR : $SGE_O_WORKDIR"
echo " nslots : $NSLOTS "
echo " TMPDIR : $TMPDIR "

cd $TMPDIR 
rsync -au $SGE_O_WORKDIR/* .

export CUDA_VISIBLE_DEVICES=`cat $TMPDIR/.gpus`



export name=BASE_NAME
#export stem="$SGE_O_WORKDIR/$name"

echo "Running Step: Equilibration 1" >> run.log
pmemd.cuda -O -i eq_1.in -p ${name}.top -c ${name}_min2.rst7 -r ${name}.eq_1.r -x ${name}.eq_1.x -e ${name}.eq_1.e -o ${name}.eq_1.o -ref ${name}_min2.rst7 -inf eq_1.inf
gzip -f ${name}.eq_1.e ${name}.eq_1.o

for i in {2..6};
do 
    echo "Running Step: Equilibration $i" >> run.log
    pmemd.cuda -O -i eq_$i.in -p ${name}.top -c ${name}.eq_$((i-1)).r -r ${name}.eq_$i.r -x ${name}.eq_$i.x -e ${name}.eq_$i.e -o ${name}.eq_$i.o -ref ${name}.eq_$((i-1)).r -inf eq_$i.inf
    gzip -f ${name}.eq_$i.e ${name}.eq_$i.o
    rsync -au * $SGE_O_WORKDIR/.
done

exit
