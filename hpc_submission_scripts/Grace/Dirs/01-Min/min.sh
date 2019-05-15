#!/bin/bash -l
#$ -N Minim
#$ -l h_rt=1:0:0 
#$ -l mem=1G
#$ -l tmpfs=10G
#$ -pe mpi 16
#$ -wd /home/YOUR USERNAME /Scratch/ DIRECTORY

module load libmatheval/1.1.11
module load plumed/2.4.3/intel-2018
module load gromacs/2018.3/plumed/intel-2018

export OMP_NUM_THREADS=1

export name=$(cd ..; basename -- "$PWD")
export GMX=/shared/ucl/apps/gromacs/2018.3/plumed/intel-2018/bin/gmx_mpi

cp ../00-Prep/EKEK.top ./$name.top
cp ../00-Prep/EKEK.gro ./$name.gro
cp ../00-Prep/EKEK_Protein*.itp . 

$GMX grompp -f min.mdp -c $name.gro -p $name.top -o min.tpr -r $name.gro

gerun $GMX mdrun -v -s min.tpr -deffnm min -tunepme -pin on -maxh 1.0

#Monitor the energy
echo 10 0 | $GMX energy -f min.edr -o ${name}_min_energy.xvg 

#produce the gro file for NVT equilibration
echo 0 | $GMX trjconv -s min.tpr -f min.trr -o eq.gro -pbc whole 


