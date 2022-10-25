#!/bin/bash -l

#SBATCH --job-name=nvt
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --constraint=gpu
#SBATCH --account=pr49

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export CRAY_CUDA_MPS=1

# Load your vesrion of Gromacs
module load daint-gpu                                                             
module use /apps/daint/UES/6.0.UP07/sandboxes/hvictor/easybuild/modules/all        
module load GROMACS/2018-CrayGNU-18.08-PLUMED-2.4.2-cuda-9.1

export FN=$(cd ..; basename -- "$PWD")
export GMX=gmx

cp ../00-Prep/posre*.itp .
cp ../01-Min/eq.gro .
cp ../01-Min/$FN.top .
cp ../00-Prep/i.ndx .
cp ../00-Prep/EKEK_Protein*.itp . 

$GMX grompp -f nvt.mdp -c eq.gro -p $FN.top -o NVT.tpr -r eq.gro -n i.ndx
srun gmx_mpi mdrun -s NVT.tpr -ntomp 1 -deffnm NVT -maxh 2

echo Total-Energy | $GMX energy -f NVT.edr -o ${FN}_NVT_energy.xvg #extract energy of system
echo Temperature | $GMX energy -f NVT.edr -o ${FN}_NVT_temperature.xvg #extract temperature of system
echo 0 | $GMX trjconv -s NVT.tpr -f NVT.trr -o NVTed.gro -b 1000 -e 1000 -pbc whole #take the last frame of NVT for the next step: NPT
echo 0 | $GMX trjconv -s NVT.tpr -f NVT.trr -o NVT_reimaged.trr -pbc whole #reimage all of  NVT simulation

