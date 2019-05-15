#$ -N Equ-NVT
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
cp ../01-Min/eq.gro .
cp ../01-Min/$FN.top .
cp ../00-Prep/EKEK_Protein*.itp . 

$GMX grompp -f nvt.mdp -c eq.gro -p $FN.top -o NVT.tpr -r eq.gro 

gerun $GMX mdrun -v -s nvt.tpr -deffnm NVT -tunepme -pin on -maxh 12.0

echo Total-Energy | $GMX energy -f NVT.edr -o ${FN}_NVT_energy.xvg #extract energy of system
echo Temperature | $GMX energy -f NVT.edr -o ${FN}_NVT_temperature.xvg #extract temperature of system
echo 0 | $GMX trjconv -s NVT.tpr -f NVT.trr -o NVTed.gro -b 1000 -e 1000 -pbc whole #take the last frame of NVT for the next step: NPT
echo 0 | $GMX trjconv -s NVT.tpr -f NVT.trr -o NVT_reimaged.trr -pbc whole #reimage all of  NVT simulation
