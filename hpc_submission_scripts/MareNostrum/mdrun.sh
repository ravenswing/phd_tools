#!/bin/bash
#SBATCH --job-name="MD_KKEE" 
#SBATCH --time=25:00:00
#SBATCH --error=MDrun_errors.%j
#SBATCH --output=MDrun_output.%j
#SBATCH --nodes=2

module load plumed/2.3.2_libmatheval
module load gromacs/5.1.4-plumed-libmatheval

export FN=$(basename -- "$PWD")
export GMX=gmx_mpi

$GMX  grompp -f mdrun.mdp -c ${FN}_1barframe.gro -p ${FN}.top -o mdrun.tpr 
mpirun $GMX  mdrun -s mdrun.tpr -deffnm MDrun -append -maxh 23.5

cp MDrun.cpt MDrun_1.cpt
sbatch mdcont.sh
