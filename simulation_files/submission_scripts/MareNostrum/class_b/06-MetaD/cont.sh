#SBATCH --output=run.out
#SBATCH --error=run.err
#SBATCH --time=24:00:00
#SBATCH --qos=class_b
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

# ------------------------------ Further MetaD --------------------------------

# set the max time
tmax=500000

# define std filenames
tpr=min.tpr
traj=metad_${FN}
ndx=i.ndx

srun $GMX mdrun -s prod.tpr -deffnm ${traj} -cpi ${traj}.cpt -append -plumed plumed_${FN}.dat -maxh 70.5

num=RUN_NUMBER 

cp ${traj}.cpt ${traj}_${num}.cpt
cp ${traj}_prev.cpt ${traj}_${num}_prev.cpt

if [ ! -f ${traj}.gro ]; then bash CONTINUE.sh;

else
    echo Protein_LIG         | $GMX trjconv -s $tpr -f ${traj}.xtc       -o ${traj}_whole.xtc -pbc whole -n i.ndx
    echo Protein Protein_LIG | $GMX trjconv -s $tpr -f ${traj}_whole.xtc -o ${traj}_clust.xtc -pbc cluster -n i.ndx
    echo Protein Protein_LIG | $GMX trjconv -s $tpr -f ${traj}_clust.xtc -o ${traj}_final.xtc -pbc mol -ur compact -center -n i.ndx
    echo Protein_LIG         | $GMX trjconv -s $tpr -f ${traj}_final.xtc -o ${FN}_short.xtc -dt 10 -n i.ndx
    echo Protein_LIG         | $GMX trjconv -s $tpr -f ${traj}_final.xtc -o ${FN}_lastframe.pdb -b $tmax -e $tmax -n i.ndx
fi
