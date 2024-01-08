#SBATCH -t 1-0
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH -p std
#SBATCH -N 1
#SBATCH -n 48 ## Edit to fit your node size; request a full node

module load apps/gromacs/2023_plumed

export GMX=gmx_mpi


# -------------------------- Start Production uMD  --------------------------

# define std filenames (before changing dir)
FN='a2b1+A769'
traj=md_${FN}
method=$(basename -- "$PWD")

$GMX grompp -f prod.mdp -c $FN.gro -p $FN.top -o prod.tpr -t md.cpt -r $FN.gro -n i.ndx -pp processed.top

## Set up the environment here according to your local configuration ($PATH, etc.)
export OMP_NUM_THREADS=4

for i in 12 8 4  ## Edit this list to fit your machine node size as per instructions (note that you'll need 4x this amount of cores)
do
    mpirun -np $i $GMX mdrun -s prod.tpr -noconfout -maxh 23.5 -g BNCHMK_$((i*4)).log
done
