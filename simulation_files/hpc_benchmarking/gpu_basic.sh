#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -N 1
#SBATCH -t 1-0
#SBATCH -e %x.%j.err
#SBATCH -o %x.%j.out

module load apps/gromacs/2023_plumed_2.9

export GMX=gmx_mpi

ncores=24


# -------------------------- Start Production uMD  --------------------------

# define std filenames (before changing dir)
FN='a2b1+A769'
traj=md_${FN}
method=$(basename -- "$PWD")

$GMX grompp -f prod.mdp -c $FN.gro -p $FN.top -o prod.tpr -t md.cpt -r $FN.gro -n i.ndx -pp processed.top

export GMX_DISABLE_GPU_TIMING=yes
OMP_NUM_THREADS=$ncores srun -n 1 -c $ncores $GMX mdrun -dlb auto -pin auto -s prod.tpr -deffnm ${traj} -maxh 23.5 -noconfout

