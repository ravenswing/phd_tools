#!/bin/bash
#SBATCH --job-name="MD2_HHHP+" 
#SBATCH --time=24:00:00
#SBATCH --error=MDrun_errors.%j
#SBATCH --output=MDrun_output.%j
#SBATCH --nodes=2

module load plumed/2.3.2_libmatheval
module load gromacs/5.1.4-plumed-libmatheval

export GMX=gmx_mpi

$GMX convert-tpr -s mdrun.tpr -extend 1000000 -o MDrun_ex.tpr
mpirun $GMX  mdrun -v -s MDrun_ex.tpr -cpi MDrun01.cpt -deffnm MDrun02 -append -maxh 23.5 

