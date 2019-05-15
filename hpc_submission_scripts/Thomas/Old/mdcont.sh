#!/bin/bash -l
#$ -l h_rt=48:0:0
#$ -l mem=1G
#$ -l tmpfs=5G 
#$ -pe mpi 48   # MUST BE A MULTIPLE OF 24
#$ -N MD_KKEE 

# !!! CHANGE TO JOB WORKING DIRECTORY !!!
#$ -wd /home/zcqsrev/baker_peptides/CHARMM36m/KKEE 

module load gcc-libs/4.9.2
module unload compilers mpi
module load compilers/intel/2017/update4
module load mpi/intel/2017/update3/intel
module load libmatheval
module load flex
module load openblas/0.2.14/intel-2015-update2
module load plumed/2.4.1/intel-2017-update4
module load gromacs/2016.4/plumed/intel-2017

export FN=$(basename -- "$PWD")
export GMX=gmx_mpi

mpirun $GMX  mdrun -s mdrun.tpr -deffnm MDrun -cpi MDrun.cpt -append -maxh 47.5

num=$((`ls -ltrh *.cpt | grep -v prev | tail -n 1 | awk '{split($NF, a, "[_.]"); print a[length(a)-1]}'`+1))

cp MDrun.cpt MDrun_${num}.cpt
cp MDrun_prev.cpt MDrun_${num}_prev.cpt

if [ ! -f MDrun.gro ]; then sbatch mdcont.sh; 

else
    echo 0 | $GMX trjconv -s mdrun.tpr -f MDrun.trr -o MDrun_lastframe.gro -b 1000000 -e 1000000 -pbc whole
    echo Protein | $GMX trjconv -s mdrun.tpr -f MDrun.trr -o ${FN}_protein.gro -b 1000000 -e 1000000 -pbc whole
    echo 0 | $GMX trjconv -s mdrun.tpr -f MDrun.trr -o MDrun_reimaged.xtc -pbc whole
    echo Protein Protein | $GMX trjconv -s min.tpr -f MDrun_reimaged.xtc -o MDrun_protraj.xtc -pbc nojump -center
    echo Protein Protein | $GMX trjconv -s min.tpr -f MDrun_protraj.xtc -o MDrun_cluster.xtc -pbc cluster
    echo Protein Protein | $GMX trjconv -s min.tpr -f MDrun_cluster.xtc -o MDrun_final.xtc -pbc mol -ur compact -center
    echo Protein | $GMX gyrate -s min.tpr -f MDrun_final.xtc -o ${FN}_Rgyr.xvg
    echo Protein | $GMX trjconv -s min.tpr -f MDrun_final.xtc -o ${FN}_protein.pdb -b 1000000 -e 1000000
    echo Backbone Backbone | $GMX rms -s min.tpr -f MDrun_final.xtc -o ${FN}_RMSD.xvg -fit rot+trans
    echo Protein | $GMX gyrate -s min.tpr -f ${FN}_initial.pdb -o ${FN}_InitialRgyr.xvg
fi
