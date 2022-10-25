#!/bin/bash -l
#$ -l h_rt=48:0:0
#$ -l mem=1G
#$ -l tmpfs=5G 
#$ -pe mpi 48   # MUST BE A MULTIPLE OF 24
#$ -N MD_MOLECULE

# !!! CHANGE TO JOB WORKING DIRECTORY !!!
#$ -wd /home/zcqsrev/Scratch/MOLECULE

module load gcc-libs/4.9.2
module unload compilers mpi
module load compilers/intel/2017/update4
module load mpi/intel/2017/update3/intel
module load libmatheval
module load flex
module load openblas/0.2.14/intel-2015-update2
module load plumed/2.4.1/intel-2017-update4
module load gromacs/2016.4/plumed/intel-2017

export FN=$(basename -- "$PWD")
export GMX=gmx_mpi

$GMX convert-tpr -s mdrun.tpr -extend 1000000 -o MDrun_ex.tpr
mpirun $GMX  mdrun -v -s MDrun_ex.tpr -cpi MDrun01.cpt -deffnm MDrun02 -append -maxh 23.5 

