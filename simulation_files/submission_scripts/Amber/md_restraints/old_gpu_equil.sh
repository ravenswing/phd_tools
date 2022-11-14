#!/bin/bash
#$ -S /bin/bash
#$ -N BASE_NAME
#$ -pe gpu 1
#$ -q fartorgpu.q
#$ -cwd 
#$ -o h.out
#$ -e e.err

# Load the modules needed 
. /etc/profile
module load amber/18_cuda9.0

# Copy inputs and files needed to the directory where the jobs will run
echo " SGE_O_WORKDIR : $SGE_O_WORKDIR"
echo " nslots : $NSLOTS "
echo " TMPDIR : $TMPDIR "
cd $TMPDIR 
rsync -au $SGE_O_WORKDIR/* .

export CUDA_VISIBLE_DEVICES=`cat $TMPDIR/.gpus`

export name=BASE_NAME

i=6

echo "Running Step: Equilibration $i" >> run.log
pmemd.cuda -O -i eq_$i.in -p ${name}.top -c ${name}.eq_$((i-1)).r -r ${name}.eq_$i.r -x ${name}.eq_$i.x -e ${name}.eq_$i.e -o ${name}.eq_$i.o -ref ${name}.eq_$((i-1)).r -inf eq_$i.inf
gzip -f ${name}.eq_$i.e ${name}.eq_$i.o

rsync -au * $SGE_O_WORKDIR/.

exit
