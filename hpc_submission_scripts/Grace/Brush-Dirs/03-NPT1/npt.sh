#$ -N Equ-NPT
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
cp ../02-NVT/NVTed.gro .
cp ../02-NVT/$FN.top .
cp ../00-Prep/EKEK_Protein*.itp . 

$GMX grompp -f npt.mdp -c NVTed.gro -p $FN.top -o NPT.tpr -r NVTed.gro

gerun $GMX mdrun -v -s NPT.tpr -deffnm NPT -tunepme -pin on -maxh 12.0

echo Total-Energy | $GMX energy -f NPT.edr -o ${FN}_NVT_energy.xvg #extract energy of system
echo Temperature | $GMX energy -f NPT.edr -o ${FN}_NVT_temperature.xvg #extract temperature of system
echo 0 | $GMX trjconv -s NPT.tpr -f NPT.trr -o NPTed.gro -b 10000 -e 10000 -pbc whole #take the last frame of NPT for the next step: MD
echo 0 | $GMX trjconv -s NPT.tpr -f NPT.trr -o NPT_reimaged.trr -pbc whole #reimage all of  NPT simulation
