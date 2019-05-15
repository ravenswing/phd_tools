#$ -N Equ-NPT2
#$ -l h_rt=12:0:0 
#$ -l mem=1G
#$ -l tmpfs=10G
#$ -pe mpi 16


module load libmatheval/1.1.11
module load plumed/2.4.3/intel-2018
module load gromacs/2018.3/plumed/intel-2018

export OMP_NUM_THREADS=1

export FN=$(cd ..; basename -- "$PWD")
export GMX=/shared/ucl/apps/gromacs/2018.3/plumed/intel-2018/bin/gmx_mpi

cp ../00-Prep/posre*.itp .
cp ../03-NPT1/NPTed.gro .
cp ../03-NPT1/$FN.top .
cp ../00-Prep/EKEK_Protein*.itp . 

head -n 14 posres_CAlpha.itp > C-alpha_pr.itp
sed -i -e '/tip4pd.itp/a \\n#ifdef CAPOSRES\n#include "Calpha_pr.itp"\n#endif\n' $FN.top

$GMX grompp -f npt2.mdp -c NPTed.gro -p $FN.top -o NPT2.tpr -r NPTed.gro

gerun $GMX mdrun -v -s NPT2.tpr -deffnm NPT2 -tunepme -pin on -maxh 12.0

echo 0 | $GMX trjconv -s NPT2.tpr -f NPT2.trr -o NPT2ed.gro -b 10000 -e 10000 -pbc whole 
echo 0 | $GMX trjconv -s NPT2.tpr -f NPT2.trr -o NPT2_reimaged.trr -pbc whole
echo 0 | $GMX trjconv -s NPT2.tpr -f NPT2.trr -o ${FN}_NPT2.xtc -pbc whole

#echo 10 0 | $GMX energy -f NPT2.edr -o pos-rest.xvg
echo 12 0 | $GMX energy -f NPT2.edr -o ${FN}_tot_energy.xvg #extract total energy of system
echo 10 0 | $GMX energy -f NPT2.edr -o ${FN}_pot_energy.xvg #extract potential energy of system
echo 11 0 | $GMX energy -f NPT2.edr -o ${FN}_kin_energy.xvg #extract kinetic energy of system
echo 13 0 | $GMX energy -f NPT2.edr -o ${FN}_temperature.xvg #extract temperature of system
echo 21 0 | $GMX energy -f NPT2.edr -o ${FN}_density.xvg #extract density of system
echo 15 0 | $GMX energy -f NPT2.edr -o ${FN}_pressure.xvg #extract pressure of system


