#$ -N Equ-NPT2
#$ -l h_rt=24:0:0 
#$ -l mem=5G
#$ -l tmpfs=10G
#$ -pe mpi 24


module load libmatheval/1.1.11
module load plumed/2.4.3/intel-2018
module load gromacs/2018.3/plumed/intel-2018

export OMP_NUM_THREADS=1

export FN=$(cd ..; basename -- "$PWD")
#export GMX=/shared/ucl/apps/gromacs/2018.3/plumed/intel-2018/bin/gmx_mpi
export GMX=gmx_mpi

cp ../00-Prep/posre*.itp .
cp ../04-NPT2/NPT2.cpt .
cp ../04-NPT2/$FN.top .
cp ../04-NPT2/NPT2ed.gro .
cp ../00-Prep/EKEK_Protein*.itp . 
cp ../01-Min/min.tpr .

head -n 7 posres_XYZAnchor.itp > XYZAnchor_pr.itp
head -n 7 posres_ZAnchor.itp > ZAnchor_pr.itp
sed -i -e '/EKEK_Protein/a \\n#ifdef ZANCHOR\n#include "ZAnchor_pr.itp"\n#endif\n ' $FN.top
sed -i -e '/EKEK_Protein/a \\n#ifdef XYZANCHOR\n#include "XYZAnchor_pr.itp"\n#endif\n ' $FN.top

gerun $GMX  grompp -f md.mdp -c NPT2ed.gro -p $FN.top -o md.tpr -t NPT2.cpt -maxwarn 1 -r NPT2ed.gro

gerun $GMX mdrun -v -s md.tpr -deffnm md -tunepme -pin on -maxh 24.0

cp md.cpt md_1.cpt
cp md_prev.cpt md_1_prev.cpt


