#!/bin/bash -l

#SBATCH --job-name=min
#SBATCH --time=00:30:00
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

export name=$(cd ..; basename -- "$PWD")
export GMX=gmx

cp ../00-Prep/topol.top ./$name.top
cp ../00-Prep/$name.gro .
cp ../00-Prep/ligand*.itp .
cp ../00-Prep/i.ndx .

$GMX grompp -f min.mdp -c $name.gro -p $name.top -o min.tpr -n i.ndx
srun gmx_mpi mdrun -s min.tpr -ntomp 1 -deffnm min -maxh 1

#Monitor the energy
echo 10 0 | $GMX energy -f min.edr -o ${name}_min_energy.xvg

#produce the gro file for NVT equilibration
echo 0 | $GMX trjconv -s min.tpr -f min.trr -o eq.gro -pbc whole

