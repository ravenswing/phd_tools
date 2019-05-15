#!/bin/bash -l

#SBATCH --job-name=npt-2
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
cp ../00-Prep/EKEK_Protein*.itp . 
cp ../00-Prep/i.ndx .
cp ../03-NPT1/NPTed.gro .
cp ../03-NPT1/$FN.top .

head -n 14 posres_CAlpha.itp > C-alpha_pr.itp   
sed -i -e '/EKEK_Protein/a \\n#ifdef CAPOSRES\n#include "C-alpha_pr.itp"\n#endif\n ' $FN.top

$GMX grompp -f npt2.mdp -c NPTed.gro -p $FN.top -o NPT2.tpr -r NPTed.gro -n i.ndx -pp processed.top
srun gmx_mpi mdrun -s NPT2.tpr -ntomp 1 -deffnm NPT2 -maxh 24


echo 0 | $GMX trjconv -s NPT2.tpr -f NPT2.trr -o NPT2ed.gro -b 10000 -e 10000 -pbc whole #take the last frame of NPT
echo 0 | $GMX trjconv -s NPT2.tpr -f NPT2.trr -o NPT2_reimaged.trr -pbc whole #reimage all of  NPT2 simulation
echo 0 | $GMX trjconv -s NPT2.tpr -f NPT2.trr -o ${FN}_NPT2.xtc -pbc whole

echo 12 0 | $GMX energy -f NPT2.edr -o ${FN}_tot_energy.xvg #extract total energy of system
echo 10 0 | $GMX energy -f NPT2.edr -o ${FN}_pot_energy.xvg #extract potential energy of system
echo 11 0 | $GMX energy -f NPT2.edr -o ${FN}_kin_energy.xvg #extract kinetic energy of system
echo 13 0 | $GMX energy -f NPT2.edr -o ${FN}_temperature.xvg #extract temperature of system
echo 21 0 | $GMX energy -f NPT2.edr -o ${FN}_density.xvg #extract density of system
echo 15 0 | $GMX energy -f NPT2.edr -o ${FN}_pressure.xvg #extract pressure of system

