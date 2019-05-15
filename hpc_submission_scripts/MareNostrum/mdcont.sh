#!/bin/bash
#SBATCH --job-name="MD_KKEE" 
#SBATCH --time=72:00:00
#SBATCH --error=MDrun_errors.%j
#SBATCH --output=MDrun_output.%j
#SBATCH --nodes=2

module load plumed/2.3.2_libmatheval
module load gromacs/5.1.4-plumed-libmatheval

export FN=$(basename -- "$PWD")
export GMX=gmx_mpi

mpirun $GMX  mdrun -s mdrun.tpr -deffnm MDrun -cpi MDrun.cpt -append -maxh 71

num=$((`ls -ltrh *.cpt | grep -v prev | tail -n 1 | awk '{split($NF, a, "[_.]"); print a[length(a)-1]}'`+1))

cp MDrun.cpt MDrun_${num}.cpt
cp MDrun_prev.cpt MDrun_${num}_prev.cpt

if [ ! -f MDrun.gro ]; then sbatch mdcont.sh; 

else
    echo 0 | $GMX trjconv -s mdrun.tpr -f MDrun.trr -o MDrun_lastframe.gro -b 1000000 -e 1000000 -pbc whole
    echo Protein | $GMX trjconv -s mdrun.tpr -f MDrun.trr -o ${FN}_protein.gro -b 1000000 -e 1000000 -pbc whole
    echo Protein | $GMX trjconv -s min.tpr -f MD_final.xtc -o ${FN}_protein.pdb -b 1000000 -e 1000000
    echo 0 | $GMX trjconv -s mdrun.tpr -f MDrun.trr -o MDrun_reimaged.xtc -pbc whole
    echo Protein Protein | $GMX trjconv -s min.tpr -f MD_reimaged.xtc -o MD_protraj.xtc -pbc nojump -center
    echo Protein Protein | $GMX trjconv -s min.tpr -f MD_protraj.xtc -o MD_cluster.xtc -pbc cluster
    echo Protein Protein | $GMX trjconv -s min.tpr -f MD_cluster.xtc -o MD_final.xtc -pbc mol -ur compact -center
    echo Protein | $GMX gyrate -s min.tpr -f MD_final.xtc -o ${FN}_Rgyr.xvg
    echo Backbone Backbone | $GMX rms -s min.tpr -f MD_final.xtc -o ${FN}_RMSD.xvg -fit rot+trans
fi
