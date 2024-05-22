#!/bin/bash
#SBATCH -J IL-8-a2b1
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -N 1
#SBATCH -t 3-0
#SBATCH -e %x.%j.err
#SBATCH -o %x.%j.out

module load apps/gromacs/2023_plumed_2.9

export GMX=gmx_mpi

# -------------------------- Start Production uMD  --------------------------

# define std filenames
FN=$(cd ..; basename -- "$PWD")
traj=metad_${FN}
method=$(cd ../..; basename -- "$PWD")
tpr=prod.tpr


$GMX grompp -f prod.mdp -c md.gro -p $FN.top -o $tpr -t md.cpt -r md.gro -n i.ndx -pp processed.top

export GMX_DISABLE_GPU_TIMING=yes

OMP_NUM_THREADS=24 srun -n 1 -c 24 gmx_mpi mdrun -dlb auto -pin auto -s prod.tpr -deffnm ${traj} -maxh 71.5 -plumed plumed_${FN}.dat

if [ ! -f ${traj}.gro ]; then sbatch _cont.sh;
fi
