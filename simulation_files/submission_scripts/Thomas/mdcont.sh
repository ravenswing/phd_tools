#$ -l h_rt=24:0:0 
#$ -l mem=5G
#$ -l tmpfs=10G
#$ -P Free
#$ -A UCL_chemM_Gervasio
#$ -pe mpi 48

module load libmatheval/1.1.11 \
	plumed/2.4.3/intel-2018 \
	gromacs/2018.3/plumed/intel-2018

export OMP_NUM_THREADS=1

export FN=$(cd ..; basename -- "$PWD")

GMX=gmx_mpi

gerun $GMX mdrun -v \
	-s MDrun_ex.tpr \
	-deffnm MDrun \
	-cpi MDrun.cpt \
	-plumed plumed \
	-maxh 23.5 \
	-append -tunepme -pin on

num=RUN_NUMBER

cp md.cpt md_${num}.cpt
cp MDrun.cpt MDrun_ext_$num.cpt
cp MDrun_prev.cpt MDrun_ext_${num}_prev.cpt


if [ ! -f MDrun.gro ]; then bash CONT.sh;

else 
    echo 0 | $GMX trjconv -s MDrun_ex.tpr -f MDrun.trr -o md_lastframe.gro -b 1500000 -e 1500000 -pbc whole
    echo Protein | $GMX trjconv -s MDrun_ex.tpr -f MDrun.trr -o ${FN}_protein.gro -b 1500000 -e 1500000 -pbc whole
    echo 0 | $GMX trjconv -s MDrun_ex.tpr -f MDrun.trr -o MD_reimaged.xtc -pbc whole
    echo Protein Protein | $GMX trjconv -s min.tpr -f MD_reimaged.xtc -o MD_protraj.xtc -pbc nojump -center
    echo Protein Protein | $GMX trjconv -s min.tpr -f MD_protraj.xtc -o MD_cluster.xtc -pbc cluster
    echo Protein Protein | $GMX trjconv -s min.tpr -f MD_cluster.xtc -o MD_final.xtc -pbc mol -ur compact -center
    echo Protein | $GMX trjconv -s min.tpr -f MD_final.xtc -o ${FN}_protein.pdb -b 1500000 -e 1500000
    echo Protein | $GMX gyrate -s min.tpr -f MD_final.xtc -o ${FN}_Rgyr.xvg
    echo Backbone Backbone | $GMX rms -s min.tpr -f MD_final.xtc -o ${FN}_RMSD.xvg -fit rot+trans
fi


