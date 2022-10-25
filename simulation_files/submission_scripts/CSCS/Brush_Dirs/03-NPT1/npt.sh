#!/bin/bash -l

#SBATCH --job-name=npt-1
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=12
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
cp ../02-NVT/NVTed.gro .
cp ../02-NVT/$FN.top .
cp ../00-Prep/i.ndx .
cp ../00-Prep/EKEK_Protein*.itp . 

$GMX grompp -f npt.mdp -c NVTed.gro -p $FN.top -o NPT.tpr -r NVTed.gro -n i.ndx
srun gmx_mpi mdrun -s NPT.tpr -ntomp 1 -deffnm NPT -maxh 24

echo Total-Energy | $GMX energy -f NPT.edr -o ${FN}_NVT_energy.xvg #extract energy of system
echo Temperature | $GMX energy -f NPT.edr -o ${FN}_NVT_temperature.xvg #extract temperature of system
echo 0 | $GMX trjconv -s NPT.tpr -f NPT.trr -o NPTed.gro -b 10000 -e 10000 -pbc whole #take the last frame of NPT 
echo 0 | $GMX trjconv -s NPT.tpr -f NPT.trr -o NPT_reimaged.trr -pbc whole #reimage all of  NPT simulation

