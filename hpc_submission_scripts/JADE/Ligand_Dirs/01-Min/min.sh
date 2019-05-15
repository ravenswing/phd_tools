#!/bin/bash 
### Each node on JADE has 8 GPUs and 40 CPUs 
### maximum walltime when running on the "big" partition (default) is 24h 
### cpus-per-task*greps = 40 
### cpus-per-task correspond to the number of CPUs assigned for each GPU task  

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --gres=gpu:1
#SBATCH --time=10:00:00 
#SBATCH -J min
#SBATCH --partition=small

module purge 
module load gromacs/2018.3 

export name=$(cd ..; basename -- "$PWD")
export GMX=gmx
export PLUMED_USE_LEPTON=yes 

cp ../00-Prep/topol.top ./$name.top # if prep_ligand used
cp ../00-Prep/complex_${name}_full.gro ./$name.gro
cp ../00-Prep/index_${name}.ndx ./i.ndx
cp ../00-Prep/lig*.itp .

$GMX grompp -f min.mdp -c $name.gro -p $name.top -o min.tpr -n i.ndx

mpirun -np ${SLURM_NTASKS_PER_NODE} --bind-to socket \
	mdrun_mpi -s min.tpr -deffnm min -maxh 10.0  \
        -ntomp ${SLURM_CPUS_PER_TASK} &> run.out

#Monitor the energy
echo 10 0 | $GMX energy -f min.edr -o ${name}_min_energy.xvg 

#produce the gro file for NVT equilibration
echo 0 | $GMX trjconv -s min.tpr -f min.trr -o eq.gro -pbc whole 


