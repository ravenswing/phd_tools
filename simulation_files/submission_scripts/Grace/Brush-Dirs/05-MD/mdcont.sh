#$ -l h_rt=24:0:0 
#$ -l mem=5G
#$ -l tmpfs=10G
#$ -pe mpi 24

module load libmatheval/1.1.11
module load plumed/2.4.3/intel-2018
module load gromacs/2018.3/plumed/intel-2018

export OMP_NUM_THREADS=1

export FN=$(cd ..; basename -- "$PWD")
export GMX=/shared/ucl/apps/gromacs/2018.3/plumed/intel-2018/bin/gmx_mpi

gerun $GMX mdrun -v -s md.tpr -deffnm md -tunepme -pin on -cpi md.cpt -append -maxh 23.5

num=RUN_NUMBER

cp md.cpt md_${num}.cpt
cp md_prev.cpt md_${num}_prev.cpt

if [ ! -f md.gro ]; then bash CONT.sh;

else 
    echo 0 | $GMX trjconv -s md.tpr -f md.trr -o md_lastframe.gro -b 500000 -e 500000 -pbc whole
    echo Protein | $GMX trjconv -s md.tpr -f md.trr -o ${FN}_protein.gro -b 500000 -e 500000 -pbc whole
    echo 0 | $GMX trjconv -s md.tpr -f md.trr -o MD_reimaged.xtc -pbc whole
    echo Protein Protein | $GMX trjconv -s min.tpr -f MD_reimaged.xtc -o MD_protraj.xtc -pbc nojump -center
    echo Protein Protein | $GMX trjconv -s min.tpr -f MD_protraj.xtc -o MD_cluster.xtc -pbc cluster
    echo Protein Protein | $GMX trjconv -s min.tpr -f MD_cluster.xtc -o MD_final.xtc -pbc mol -ur compact -center
    echo Protein | $GMX trjconv -s min.tpr -f MD_final.xtc -o ${FN}_protein.pdb -b 500000 -e 500000
    echo Protein | $GMX gyrate -s min.tpr -f MD_final.xtc -o ${FN}_Rgyr.xvg
    echo Backbone Backbone | $GMX rms -s min.tpr -f MD_final.xtc -o ${FN}_RMSD.xvg -fit rot+trans
fi


