#!/bin/bash
### Each node on JADE has 8 GPUs and 40 CPUs 

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=5
#SBATCH --gres=gpu:8
#SBATCH --time=24:00:00 
#SBATCH -J npt-5alo
#SBATCH --partition=big

module purge 
module load gromacs/2018.3 

export FN=$(cd ..; basename -- "$PWD")
export GMX=gmx

cp ../00-Prep/posre*.itp .
cp ../00-Prep/index_$FN.ndx ./i.ndx
cp ../01-Min/lig*.itp .
cp ../02-NVT/NVTed.gro .
cp ../02-NVT/$FN.top .

$GMX grompp -f npt.mdp -c NVTed.gro -p $FN.top -o NPT.tpr -r NVTed.gro -n i.ndx

#mpirun -np ${SLURM_NTASKS_PER_NODE} --bind-to socket \
#	mdrun_mpi -s NPT.tpr -deffnm NPT -maxh 24.0  \
#        -ntomp ${SLURM_CPUS_PER_TASK} &> run.out

mpirun -np ${SLURM_NTASKS_PER_NODE} --bind-to socket mdrun_mpi -s NPT.tpr -v -pin on -nb gpu -tunepme -pme gpu -pmefft gpu -npme 1 -deffnm NPT -ntomp ${SLURM_CPUS_PER_TASK} -maxh 23.5 &> run.out

echo Total-Energy | $GMX energy -f NPT.edr -o ${FN}_NVT_energy.xvg #extract energy of system
echo Temperature | $GMX energy -f NPT.edr -o ${FN}_NVT_temperature.xvg #extract temperature of system
echo 0 | $GMX trjconv -s NPT.tpr -f NPT.trr -o NPTed.gro -b 10000 -e 10000 -pbc whole #take the last frame of NPT for the next step: MD
echo 0 | $GMX trjconv -s NPT.tpr -f NPT.trr -o NPT_reimaged.trr -pbc whole #reimage all of  NPT simulation
