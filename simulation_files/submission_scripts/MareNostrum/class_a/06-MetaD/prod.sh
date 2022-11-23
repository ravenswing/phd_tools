#SBATCH --output=run.out
#SBATCH --error=run.err
#SBATCH --time=72:00:00
#SBATCH --qos=class_a
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-node=24
#SBATCH --nodes=4

module purge
module load intel/2020.1
module load impi/2018.4
module load mkl/2020.1
module load boost/1.75.0
module load plumed/2.8.0
module load gromacs/2021.4-plumed.2.8.0

export FN=$(cd ..; basename -- "$PWD")
export GMX=gmx_mpi

cp ../../../${FN}/0345-EQ-MD/md.gro .
cp ../../../${FN}/0345-EQ-MD/md.cpt .
cp ../../../${FN}/0345-EQ-MD/${FN}.top .

cp ../../../${FN}/01-Min/i.ndx .
cp ../../../${FN}/01-Min/*.itp .
cp ../../../${FN}/01-Min/min.tpr .

#sed -i 's/Water_and_ions/Water/g' prod.mdp

# -------------------------- Start Producion MetaD  --------------------------

# define std filenames
traj=metad_${FN}

$GMX grompp -f prod.mdp -c md.gro -p $FN.top -o prod.tpr -t md.cpt -r md.gro -n i.ndx -pp processed.top
srun $GMX mdrun -s prod.tpr -deffnm ${traj} -plumed plumed_${FN}.dat -maxh 71.5

cp ${traj}.cpt ${traj}_1.cpt
cp ${traj}_prev.cpt ${traj}_1_prev.cpt

