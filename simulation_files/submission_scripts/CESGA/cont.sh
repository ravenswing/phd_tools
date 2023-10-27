#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -t 3-0
#SBATCH --gres=gpu:a100:1   # Request 1 GPU of 2 available on an average A100 node
#SBATCH -c 32               # Cores per task requested
#SBATCH --mem-per-cpu=3G

module load cesga/2020 gcc/system openmpi/4.1.4_ft3_cuda gromacs/2021.4-plumed-2.8.0

export GMX=gmx_mpi


# -------------------------- Start Producion MetaD  --------------------------

# set the max time
tmax=500000
# define std filenames
FN=$(cd ..; basename -- "$PWD")
traj=metad_${FN}
method=$(cd ../..; basename -- "$PWD")
tpr=prod.tpr
ndx=i.ndx

srun $GMX mdrun -s prod.tpr -deffnm ${traj} -cpi ${traj}.cpt -append -plumed ${FN}_${method}.dat -maxh 71.5

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
