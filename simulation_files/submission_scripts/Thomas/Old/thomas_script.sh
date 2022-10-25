#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=0:05:0
#$ -l mem=1G
#$ -l tmpfs=5G
#$ -N bench
#$ -pe mpi 192

# !!! CHANGE TO JOB WORKING DIRECTORY !!!
#$ -wd /home/uccagm1/Scratch/benchmark

module load gcc-libs/4.9.2
module load compilers/intel/2017/update1
module load mpi/intel/2017/update1/intel
module load libmatheval/1.1.11
module load flex/2.5.39
module load openblas/0.2.14/intel-2015-update2
module load plumed/2.3.1/intel-2017-update1
module load gromacs/2016.3/plumed/intel-2017-update1



gmx_mpi grompp -c system.gro -f md.mdp -p system.top -o topol.tpr -n index.ndx -maxwarn 1
mpirun gmx_mpi mdrun -s topol.tpr -deffnm test -nsteps 50000 -tunepme -plumed plumed.dat

