#!/bin/bash
#SBATCH --job-name="MDrun_4000"
#SBATCH --time=25:00:00
#SBATCH --error=MDrun_errors.%j
#SBATCH --output=MDrun_output.%j
#SBATCH --nodes=8

module load plumed/2.3.2_libmatheval
module load gromacs/5.1.4-plumed-libmatheval

export GMX=gmx_mpi

$GMX  grompp -f mdrun.mdp -c 4000ext_8A_1bar.gro -p 4000_topol.top -o mdrun.tpr -maxwarn 1
mpirun $GMX  mdrun -s mdrun.tpr -deffnm MDrun

echo 0 | $GMX trjconv -s mdrun.tpr -f mdrun.trr -o MDrun.gro -b 1000 -e 1000 -pbc whole
echo 1 | $GMX trjconv -s mdrun.tpr -f mdrun.trr -o MDrun_protein.pdb -b 1000 -e 1000 -pbc whole
echo 0 | $GMX trjconv -s mdrun.tpr -f mdrun.trr -o MDrun_reimaged.xtc -pbc whole
