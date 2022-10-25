#SBATCH --time=24:00:00
#SBATCH --nodes=4
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

cp ../02-NVT/i.ndx .
cp ../01-Min/min.tpr .
cp ../04-NPT2/NPT2.cpt .
cp ../04-NPT2/$FN.top .
cp ../04-NPT2/NPT2ed.gro .
cp ../00-Prep/posre*.itp .
cp ../00-Prep/EKEK_Protein*.itp . 

head -n 7 posres_XYZAnchor.itp > XYZAnchor_pr.itp   
head -n 7 posres_ZAnchor.itp > ZAnchor_pr.itp   
sed -i -e '/EKEK_Protein/a \\n#ifdef ZANCHOR\n#include "ZAnchor_pr.itp"\n#endif\n ' $FN.top
sed -i -e '/EKEK_Protein/a \\n#ifdef XYZANCHOR\n#include "XYZAnchor_pr.itp"\n#endif\n ' $FN.top

$GMX  grompp -f md.mdp -c NPT2ed.gro -p $FN.top -o md.tpr -t NPT2.cpt -r NPT2ed.gro -n i.ndx -pp processed.top
srun gmx_mpi mdrun -s md.tpr -ntomp 1 -deffnm md -maxh 24

cp md.cpt md_1.cpt
cp md_prev.cpt md_1_prev.cpt

