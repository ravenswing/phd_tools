#!/bin/csh 
#$ -S /bin/bash
#$ -N 
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
cp ../02-NVT/NVT.gro .
cp ../02-NVT/$FN.top .
cp ../01-Min/i.ndx .
cp ../01-Min/*.itp .

# Get file names from directory structure 
FN=$(cd ..; basename -- "$PWD")
ndx=i.ndx

# Copy inputs and files needed to the directory where the jobs will run
cd $TMPDIR 
rsync -au $SGE_O_WORKDIR/* .
# Change inputs if running without ion conc. and neutral system
#sed -i 's/Water_and_ions/Water/g' npt.mdp
#sed -i 's/Water_and_ions/Water/g' npt2.mdp
#sed -i 's/Water_and_ions/Water/g' md.mdp

# Change inputs for an apo system
sed -i 's/Protein_LIG/Protein/g' npt.mdp
sed -i 's/Protein_LIG/Protein/g' npt2.mdp
sed -i 's/Protein_LIG/Protein/g' md.mdp

# ------------------------------ NVT: Berendsen ------------------------------

$GMX grompp -f npt.mdp -c NVT.gro -p $FN.top -o NPT.tpr -r NVT.gro -n $ndx
$GMX mdrun -s NPT.tpr -deffnm NPT
rsync -au * $SGE_O_WORKDIR/.

echo Total-Energy | $GMX energy -f NPT.edr -o ${FN}_NPT_energy.xvg      # extract energy of system
echo Temperature  | $GMX energy -f NPT.edr -o ${FN}_NPT_temperature.xvg # extract temperature of system
rsync -au * $SGE_O_WORKDIR/.

# -------------------------- NPT2: Parinello-Rahman --------------------------

$GMX grompp -f npt2.mdp -c NPT.gro -p $FN.top -o NPT2.tpr -r NPT.gro -n $ndx
$GMX mdrun -s NPT2.tpr -deffnm NPT2
rsync -au * $SGE_O_WORKDIR/.

echo 12 0 | $GMX energy -f NPT2.edr -o ${FN}_tot_energy.xvg  		# extract total energy of system
echo 13 0 | $GMX energy -f NPT2.edr -o ${FN}_temperature.xvg 		# extract temperature of system
echo 15 0 | $GMX energy -f NPT2.edr -o ${FN}_pressure.xvg    		# extract pressure of system
rsync -au * $SGE_O_WORKDIR/.


# -------------------------------- Initial MD --------------------------------

$GMX grompp -f md.mdp -c NPT2.gro -p $FN.top -o md.tpr -t NPT2.cpt -r NPT2.gro -n $ndx 
$GMX mdrun -s md.tpr -deffnm md -maxh 20
rsync -au * $SGE_O_WORKDIR/.

cp md.cpt md_1.cpt
cp md_prev.cpt md_1_prev.cpt
rsync -au * $SGE_O_WORKDIR/.
cd $SGE_O_WORKDIR
touch *
exit
