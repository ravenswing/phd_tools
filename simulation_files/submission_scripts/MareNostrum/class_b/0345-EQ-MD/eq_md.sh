#SBATCH --output=run.out
#SBATCH --error=run.err
#SBATCH --time=24:00:00
#SBATCH --qos=class_b
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

export FN=$(cd ..; basename -- "$PWD")
export GMX=gmx_mpi

cp ../02-NVT/NVT.gro .
cp ../02-NVT/$FN.top .
cp ../01-Min/i.ndx .
cp ../01-Min/*.itp .
cp ../01-Min/min.tpr .

#sed -i 's/Water_and_ions/Water/g' npt.mdp
#sed -i 's/Water_and_ions/Water/g' npt2.mdp
#sed -i 's/Water_and_ions/Water/g' md.mdp

# ------------------------------ NVT: Berendsen ------------------------------

$GMX grompp -f npt.mdp -c NVT.gro -p $FN.top -o NPT.tpr -r NVT.gro -n i.ndx
srun gmx_mpi mdrun -s NPT.tpr -deffnm NPT

echo Total-Energy | $GMX energy -f NPT.edr -o ${FN}_NPT_energy.xvg      #extract energy of system
echo Temperature  | $GMX energy -f NPT.edr -o ${FN}_NPT_temperature.xvg #extract temperature of system


# -------------------------- NPT2: Parinello-Rahman --------------------------

$GMX grompp -f npt2.mdp -c NPT.gro -p $FN.top -o NPT2.tpr -r NPT.gro -n i.ndx
srun gmx_mpi mdrun -s NPT2.tpr -deffnm NPT2

echo 12 0 | $GMX energy -f NPT2.edr -o ${FN}_tot_energy.xvg  #extract total energy of system
echo 13 0 | $GMX energy -f NPT2.edr -o ${FN}_temperature.xvg #extract temperature of system
echo 15 0 | $GMX energy -f NPT2.edr -o ${FN}_pressure.xvg    #extract pressure of system


# -------------------------------- Initial MD --------------------------------

$GMX  grompp -f md.mdp -c NPT2.gro -p $FN.top -o md.tpr -t NPT2.cpt -r NPT2.gro -n i.ndx -pp processed.top
srun gmx_mpi mdrun -s md.tpr -deffnm md -maxh 71.5

cp md.cpt md_1.cpt
cp md_prev.cpt md_1_prev.cpt

