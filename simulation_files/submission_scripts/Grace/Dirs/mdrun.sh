#!/bin/bash -l
#$ -N MDrun 
#$ -l h_rt=24:0:0 
#$ -l mem=1G
#$ -l tmpfs=10G
#$ -pe mpi 64 
#$ -wd /home/zcqsrev/Scratch/HHHG

export OMP_NUM_THREADS=1

module unload compilers mpi
module load compilers/intel/2015/update2
module load mpi/intel/2015/update3/intel

module load libmatheval
module load flex
module load openblas/0.2.14/intel-2015-update2
module load plumed/2.2.3/intel-2015-update2
module load gromacs/5.1.3/plumed/intel-2015-update2

export FN=$(basename -- "$PWD")
export GMX=/shared/ucl/apps/gromacs/5.1.3/plumed/intel-2015-update2/bin/gmx_mpi

$GMX  grompp -f mdrun.mdp -c ${FN}_1barframe.gro -p ${FN}.top -o mdrun.tpr -maxwarn 1
gerun $GMX  mdrun -v -s mdrun.tpr -deffnm MDrun -tunepme -pin on -cpi -append -maxh 23.5 # -npme
#gerun $GROMACS mdrun -v -s md_0_1.tpr 

echo 0 | $GMX trjconv -s mdrun.tpr -f mdrun.trr -o ${FN}_MDrun.gro -b 1000 -e 1000 -pbc whole
echo 1 | $GMX trjconv -s mdrun.tpr -f mdrun.trr -o ${FN}_MDrun_protein.pdb -b 1000 -e 1000 -pbc whole
echo 0 | $GMX trjconv -s mdrun.tpr -f mdrun.trr -o ${FN}_MDrun_reimaged.xtc -pbc whole
