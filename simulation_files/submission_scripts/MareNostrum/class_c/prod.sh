#SBATCH --output=run.out
#SBATCH --error=run.err
#SBATCH --time=24:00:00
#SBATCH --qos=class_c
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-node=24
#SBATCH --nodes=2

module purge
module load intel/2020.1
module load impi/2018.4
module load mkl/2020.1
module load boost/1.75.0
module load plumed/2.8.0
module load gromacs/2021.4-plumed.2.8.0

export GMX=gmx_mpi

# -------------------------- Start Producion MetaD  --------------------------

# define std filenames
FN=$(basename -- "$PWD")
traj=metad_${FN}
method=$(cd ../..; basename -- "$PWD")

$GMX grompp -f prod.mdp -c $FN.gro -p $FN.top -o prod.tpr -t md.cpt -r $FN.gro -n i.ndx -pp processed.top
srun $GMX mdrun -s prod.tpr -deffnm ${traj} -plumed ${FN}_${method}.dat -maxh 23.5

cp ${traj}.cpt ${traj}_1.cpt
cp ${traj}_prev.cpt ${traj}_1_prev.cpt

