#SBATCH --output=run.out
#SBATCH --error=run.err
#SBATCH --time=72:00:00
#SBATCH --qos=class_a
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-node=24
#SBATCH --nodes=4

module purge
module load intel/2020.1
module load impi/2018.4
module load mkl/2020.1
module load boost/1.75.0
module load plumed/2.8.0
module load gromacs/2021.4-plumed.2.8.0

export FN=$(cd ..; basename -- "$PWD")
export GMX=gmx_mpi

# -------------------------- Start Producion MetaD  --------------------------

srun gmx_mpi mdrun -s prod2.tpr -deffnm metad_${FN} -cpi metad_${FN}.cpt -append -plumed plumed_${FN}.dat -maxh 70.5

tpr=min.tpr
traj=metad_${FN}
ndx=i.ndx

num=RUN_NUMBER 

cp ${traj}.cpt ${traj}_${num}.cpt
cp ${traj}_prev.cpt ${traj}_${num}_prev.cpt

if [ ! -f ${traj}.gro ]; then bash CONTINUE.sh;

else
	echo Protein_LIG         | $GMX trjconv -s $tpr -f ${traj}.xtc       -o ${traj}_whole.xtc -pbc whole -n i.ndx
	echo Protein Protein_LIG | $GMX trjconv -s $tpr -f ${traj}_whole.xtc -o ${traj}_clust.xtc -pbc cluster -n i.ndx
	echo Protein Protein_LIG | $GMX trjconv -s $tpr -f ${traj}_clust.xtc -o ${traj}_final.xtc -pbc mol -ur compact -center -n i.ndx
	echo Protein_LIG         | $GMX trjconv -s $tpr -f ${traj}_final.xtc -o ${FN}_short.xtc -dt 10 -n i.ndx
	echo Protein_LIG         | $GMX trjconv -s $tpr -f ${traj}_final.xtc -o ${FN}_lastframe.pdb -b 500000 -e 500000 -n i.ndx
fi




