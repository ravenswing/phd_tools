#!/bin/bash
### Each node on JADE has 8 GPUs and 40 CPUs 
### maximum walltime when running on the "big" partition (default) is 24h 
### cpus-per-task*greps = 40 
### cpus-per-task correspond to the number of CPUs assigned for each GPU task  

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00 
#SBATCH -J npt2
#SBATCH --partition=small

module purge 
module load gromacs/2018.3 

export FN=$(cd ..; basename -- "$PWD")
export GMX=gmx

cp ../00-Prep/posre*.itp .
cp ../00-Prep/index_$FN.ndx ./i.ndx
cp ../01-Min/lig*.itp .
cp ../03-NPT1/NPTed.gro .
cp ../03-NPT1/$FN.top .
sed -i 's/posre.itp/posres_CAlpha.itp/g' $FN.top # change to position restraint selection

echo C-alpha | $GMX genrestr -f NPTed.gro -o posres_CAlpha.itp

$GMX grompp -f npt2.mdp -c NPTed.gro -p $FN.top -o NPT2.tpr -r NPTed.gro -n i.ndx

mpirun -np ${SLURM_NTASKS_PER_NODE} --bind-to socket \
	mdrun_mpi -s NPT2.tpr -deffnm NPT2 -maxh 12.0  \
        -ntomp ${SLURM_CPUS_PER_TASK} &> run.out

echo 0 | $GMX trjconv -s NPT2.tpr -f NPT2.trr -o NPT2ed.gro -b 50 -e 50 -pbc whole #take the last frame of NPT2 for the next step: MD
echo 0 | $GMX trjconv -s NPT2.tpr -f NPT2.trr -o NPT2_reimaged.trr -pbc whole #reimage all of  NPT2 simulation

echo 0 | $GMX trjconv -s NPT2.tpr -f NPT2.trr -o ${FN}_NPT2.xtc -pbc whole


#echo 10 0 | $GMX energy -f NPT2.edr -o pos-rest.xvg
echo 12 0 | $GMX energy -f NPT2.edr -o ${FN}_tot_energy.xvg #extract total energy of system
echo 10 0 | $GMX energy -f NPT2.edr -o ${FN}_pot_energy.xvg #extract potential energy of system
echo 11 0 | $GMX energy -f NPT2.edr -o ${FN}_kin_energy.xvg #extract kinetic energy of system
echo 13 0 | $GMX energy -f NPT2.edr -o ${FN}_temperature.xvg #extract temperature of system
echo 21 0 | $GMX energy -f NPT2.edr -o ${FN}_density.xvg #extract density of system
echo 15 0 | $GMX energy -f NPT2.edr -o ${FN}_pressure.xvg #extract pressure of system


