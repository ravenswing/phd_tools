#!/bin/bash
#SBATCH --chdir=/home/ub183/ub183944/scratch/ampk_replicas/a2b1+A769/0345-EQ-MD
#SBATCH --job-name=pp
#SBATCH --output=postp.out
#SBATCH --error=postp.err
#SBATCH --time=24:00:00
#SBATCH --qos=class_a
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-node=24
#SBATCH --nodes=1

module purge
module load intel/2020.1
module load impi/2018.4
module load mkl/2020.1
module load boost/1.75.0
module load plumed/2.8.0
module load gromacs/2021.4-plumed.2.8.0

export FN=$(cd ..; basename -- "$PWD")
export GMX=gmx_mpi

tpr=min.tpr
traj=md
ndx=i.ndx

echo Protein_LIG         | $GMX trjconv -s $tpr -f ${traj}.xtc       -o ${traj}_whole.xtc -pbc whole -n i.ndx
echo Protein Protein_LIG | $GMX trjconv -s $tpr -f ${traj}_whole.xtc -o ${traj}_clust.xtc -pbc cluster -n i.ndx
echo Protein Protein_LIG | $GMX trjconv -s $tpr -f ${traj}_clust.xtc -o ${traj}_final.xtc -pbc mol -ur compact -center -n i.ndx
echo Protein_LIG         | $GMX trjconv -s $tpr -f ${traj}_final.xtc -o ${FN}_lastframe.pdb -b 10000 -e 10000 -n i.ndx
