#SBATCH --job-name=nvt
#SBATCH --output=nvt.out
#SBATCH --error=nvt.err
#SBATCH --time=02:00:00
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

cp ../01-Min/min_out.gro .
cp ../01-Min/$FN.top .
cp ../01-Min/*.itp . 
cp ../01-Min/i.ndx .

#sed -i 's/Water_and_ions/Water/g' nvt.mdp

$GMX grompp -f nvt.mdp -c min_out.gro -p $FN.top -o NVT.tpr -r min_out.gro -n i.ndx
srun gmx_mpi mdrun -s NVT.tpr -deffnm NVT -maxh 2

#extract energy of system
echo Total-Energy | $GMX energy -f NVT.edr -o ${FN}_NVT_energy.xvg 

#extract temperature of system
echo Temperature | $GMX energy -f NVT.edr -o ${FN}_NVT_temperature.xvg 

#take the last frame of NVT for the next step: NPT
echo 0 | $GMX trjconv -s NVT.tpr -f NVT.trr -o NVTed.gro -b 1000 -e 1000 -pbc whole

