#!/bin/csh 
#$ -S /bin/bash
#$ -N 02-NVT
#$ -pe smp 1
#$ -q iqtc10_g1819_gpu.q
#$ -cwd 
#$ -o $JOB_NAME.out
#$ -e $JOB_NAME.err

# --------------------------------- SETUP  ------------------------------------

# Load the modules and set environment variables
source /etc/profile
module load gromacs/2021.4_cuda
export OMP_NUM_THREADS=1
export SGE_ROOT="/sge"
export GMX=gmx_mpi
export CUDA_VISIBLE_DEVICES=`cat $TMPDIR/.gpus`

# Put run information to start log file
echo " SGE_O_WORKDIR : $SGE_O_WORKDIR" > run.log
echo " nslots : $NSLOTS " >> run.log
echo " TMPDIR : $TMPDIR " >> run.log

# Copy from previous step
cp ../01-Min/min.gro .
cp ../01-Min/$FN.top .
cp ../01-Min/*.itp . 
cp ../01-Min/i.ndx .

# Get file names from directory structure
FN=$(cd ..; basename -- "$PWD")
ndx=i.ndx

# Copy inputs and files needed to the directory where the jobs will run
cd $TMPDIR 
rsync -au $SGE_O_WORKDIR/* .

#sed -i 's/Water_and_ions/Water/g' nvt.mdp

$GMX grompp -f nvt.mdp -c min.gro -p $FN.top -o NVT.tpr -r min.gro -n $ndx
$GMX mdrun -s NVT.tpr -deffnm NVT
rsync -au * $SGE_O_WORKDIR/.

# Extract energy of system
echo Total-Energy | $GMX energy -f NVT.edr -o ${FN}_NVT_energy.xvg 

# Extract temperature of system
echo Temperature | $GMX energy -f NVT.edr -o ${FN}_NVT_temperature.xvg 

# Take the last frame of NVT for visualisation
echo 0 | $GMX trjconv -s NVT.tpr -f NVT.trr -o NVTed.gro -b 1000 -e 1000 -pbc whole

# Final file copy
rsync -au * $SGE_O_WORKDIR/.
cd $SGE_O_WORKDIR
touch *
exit
