#!/bin/csh 
#$ -S /bin/bash
#$ -N BASE_NAME
#$ -pe smp 1
#$ -q iqtc10_g1819_gpu.q
#$ -cwd 
#$ -o BASE_NAME_md.out
#$ -e BASE_NAME_md.err

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

pmemd.cuda -O -i md.in -p ${name}.top -c ${name}_eq_6.r -r ${name}_md_1.r -x ${name}_md_1.x -e ${name}_md_1.e -o ${name}_md_1.o -ref ${name}_eq_6.r -inf md_1.inf
gzip -f ${name}_md_1.e ${name}_md_1.o
rsync -au * $SGE_O_WORKDIR/.

for i in {2..200..1}
do
	echo "Running Step: $i" >> run.log
	pmemd.cuda -O -i md.in -p ${name}.top -c ${name}_md_$((i-1)).r -r ${name}_md_$i.r -x ${name}_md_$i.x -e ${name}_md_$i.e -o ${name}_md_$i.o -ref  ${name}_md_$((i-1)).r -inf md_$i.inf
	gzip -f ${name}_md_$i.e ${name}_md_$i.o
	rsync -au * $SGE_O_WORKDIR/.
done

exit
