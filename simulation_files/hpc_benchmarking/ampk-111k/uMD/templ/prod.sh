#SBATCH -o %j_%i.out
#SBATCH -e %j_%i.err
#SBATCH -t 02:00:00
#SBATCH --gres=gpu:a100:1   # Request 1 GPU of 2 available on an average A100 node
#SBATCH -c 32               # Cores per task requested
#SBATCH --mem-per-cpu=3G

module load cesga/2020 gcc/system openmpi/4.1.4_ft3_cuda gromacs/2021.4-plumed-2.8.0

export GMX=gmx_mpi


# -------------------------- Start Production uMD  --------------------------

# define std filenames
FN='a2b1+A769'
traj=md_${FN}
method=$(basename -- "$PWD")

$GMX grompp -f prod.mdp -c $FN.gro -p $FN.top -o prod.tpr -t md.cpt -r $FN.gro -n i.ndx -pp processed.top
srun $GMX mdrun -s prod.tpr -deffnm ${traj} -maxh 23.5
