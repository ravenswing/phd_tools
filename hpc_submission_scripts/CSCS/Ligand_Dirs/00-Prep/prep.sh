#!/bin/bash -l

#SBATCH --job-name=prep
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --partition=prepost
#SBATCH --constraint=mc
#SBATCH --account=pr49

#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export CRAY_CUDA_MPS=1

# Load your vesrion of Gromacs
module load daint-gpu                                                             
module use /apps/daint/UES/6.0.UP07/sandboxes/hvictor/easybuild/modules/all        
module load GROMACS/2018-CrayGNU-18.08-PLUMED-2.4.2-cuda-9.1

export name=$(cd ..; basename -- "$PWD")
export GMX=gmx

# ---------- MANUAL MERGING ----------

# generate index file
echo -e " 1 | 13 \n q" | gmx make_ndx -f ${name}_edited.gro -o i.ndx

# build dodecahedral box w/ 1.2 nm clearance
$GMX editconf -f ${name}_edited.gro -o ${name}_box.gro -c -d 1.0 -bt dodecahedron

$GMX solvate -cp ${name}_box.gro -cs tip4p.gro -o ${name}_sol.gro -p topol.top          # TIP4P / TIP4P-D
##$GMX solvate -cp ${name}_box.gro -cs spc216.gro -o ${name}_sol.gro -p topol.top       # TIP3P

# USE: amber-disp       NOT: amber14sb
sed -i -e '/SOL/ s/HW2/HW3/ ' ${name}_sol.gro
sed -i -e '/SOL/ s/HW1/HW2/ ' ${name}_sol.gro
sed -i -e '/SOL/ s/ MW/MW4/ ' ${name}_sol.gro

# add ions to reach 0.15 moldm-3 NaCl concentration
$GMX grompp -f prep.mdp -c ${name}_sol.gro -p topol.top -o ions.tpr -n i.ndx -maxwarn 1
echo SOL | $GMX genion -s ions.tpr -o ${name}.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

# reshape to make viewable in 
$GMX grompp -f prep.mdp -c $name.gro -p topol.top -o reshape.tpr -n i.ndx
echo System | $GMX trjconv -f $name.gro -s reshape.tpr -o Readable${name}.gro -pbc mol -ur compact

# generate Alpha Carbon position restraints
echo C-alpha | $GMX genrestr -f $name.gro -o posres_CAlpha.itp
echo 0 | $GMX genrestr -f ${name}_ligand_crystal.gro -o posre_lig.itp
echo Protein | $GMX genrestr -f $name.gro -o posre_prot.itp

# generate index file
echo -e " 1 | 13 \n q" | gmx make_ndx -f ${name}.gro -o i.ndx

