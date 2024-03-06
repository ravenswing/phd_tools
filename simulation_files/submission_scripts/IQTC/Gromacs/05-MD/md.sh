#SBATCH --output=run.out
#SBATCH --error=run.err
#SBATCH --time=72:00:00
#SBATCH --qos=class_a
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

export FN=$(basename -- "$PWD")
export GMX=gmx_mpi


# -------------------------------- Initial MD --------------------------------

$GMX  grompp -f md.mdp -c $FN.gro -p $FN.top -o md.tpr -r $FN.gro -n i.ndx -pp processed.top
srun gmx_mpi mdrun -s md.tpr -deffnm md -maxh 71.5

cp md.cpt md_1.cpt
cp md_prev.cpt md_1_prev.cpt

