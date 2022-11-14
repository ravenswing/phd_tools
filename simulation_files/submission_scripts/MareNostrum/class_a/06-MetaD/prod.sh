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

cp ../0345-EQ-MD/md.gro .
cp ../0345-EQ-MD/md.cpt .
cp ../0345-EQ-MD/${FN}.top .

cp ../01-Min/i.ndx .
cp ../01-Min/*.itp .
cp ../01-Min/min.tpr .

#sed -i 's/Water_and_ions/Water/g' prod.mdp

# -------------------------- Start Producion MetaD  --------------------------

$GMX grompp -f prod.mdp -c md.gro -p $FN.top -o prod.tpr -t md.cpt -r md.gro -n i.ndx -pp processed.top
srun gmx_mpi mdrun -s prod.tpr -deffnm metad_${FN} -plumed plumed_${FN}.dat -maxh 71.5

cp metad_${FN}.cpt metad_${FN}_1.cpt
cp metad_${FN}_prev.cpt metad_${FN}_1_prev.cpt

