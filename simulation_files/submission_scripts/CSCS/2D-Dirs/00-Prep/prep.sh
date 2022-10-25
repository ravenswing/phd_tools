#!/bin/bash -l
#
# It is required to specify either --constraint=gpu or --constraint=mc: the option --constraint=gpu makes sure that the scheduler allocates 
# the XC50 Intel Haswell 12-core nodes with GPU devices and automatically sets the option --gres=gpu:1. 
# The --constraint=mc would ensure instead that the batch job is allocated to the multicore XC40 Intel Broadwell 2 x 18-core nodes and not to the GPU nodes.
# please note that you need to use the latter with the partition prepost and if you target large memory nodes using the flag --mem=120GB.
#
#SBATCH --job-name=lin-para
#SBATCH --time=01:00:00
#SBATCH --nodes=1
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

export name=$(cd ..; basename -- "$PWD")
export GMX=gmx

for i in {1..2}
do 
    num=$(($i * 16))
    echo $num
    sed -i -e "/HE2 HIE    $num/a TER" ${name}_initial.pdb
done


# convert pdb
echo 1 | $GMX pdb2gmx -f ${name}_initial.pdb -o ${name}_built.gro -p $name.top -ff amber_disp -ignh

# build dodecahedral box w/ 0.8 nm clearance
$GMX editconf -f ${name}_built.gro -o ${name}_box.gro -c -d 1.2 -bt dodecahedron

# solvate with TIP4P-D water
$GMX solvate -cp ${name}_box.gro -cs tip4p.gro -o ${name}_sol.gro -p $name.top

sed -i -e '/SOL/ s/HW2/HW3/ ' ${name}_sol.gro
sed -i -e '/SOL/ s/HW1/HW2/ ' ${name}_sol.gro
sed -i -e '/SOL/ s/ MW/MW4/ ' ${name}_sol.gro

$GMX grompp -f prep.mdp -c ${name}_sol.gro -p $name.top -o ions.tpr
# neutralise
#echo SOL | $GMX genion -s ions.tpr -o ${name}.gro -p $name.top -pname NA -nname CL -neutral

# add ions to reach 0.15 moldm-3 NaCl concentration
echo SOL | $GMX genion -s ions.tpr -o ${name}.gro -p $name.top -pname NA -nname CL -neutral -conc 0.15

# reshape to make viewable in 
$GMX grompp -f prep.mdp -c $name.gro -p $name.top -o reshape.tpr
echo System | $GMX trjconv -f $name.gro -s reshape.tpr -o Readable${name}.gro -pbc mol -ur compact

# generate Alpha Carbon position restraints
echo C-alpha | $GMX genrestr -f $name.gro -o posres_CAlpha.itp


