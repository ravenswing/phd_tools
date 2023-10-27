#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -t 3-0
#SBATCH --gres=gpu:a100:1   # Request 1 GPU of 2 available on an average A100 node
#SBATCH -c 32               # Cores per task requested
#SBATCH --mem-per-cpu=3G

module load cesga/2020 gcc/system openmpi/4.1.4_ft3_cuda gromacs/2021.4-plumed-2.8.0

export GMX=gmx_mpi


# -------------------------- Start Producion MetaD  --------------------------

# define std filenames
FN=$(cd ..; basename -- "$PWD")
traj=metad_${FN}
method=$(cd ../..; basename -- "$PWD")

$GMX grompp -f prod.mdp -c md.gro -p $FN.top -o prod.tpr -t md.cpt -r md.gro -n i.ndx -pp processed.top
srun $GMX mdrun -s prod.tpr -deffnm ${traj} -plumed ${FN}_${method}.dat -maxh 71.5
