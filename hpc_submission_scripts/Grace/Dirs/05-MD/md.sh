#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --gres=gpu:1
#SBATCH --time=72:00:00 
#SBATCH --partition=small

module purge 
module load gromacs/2018.3 

export FN=$(cd ..; basename -- "$PWD")
export GMX=gmx

cp ../04-NPT2/NPT2.cpt .
cp ../04-NPT2/$FN.top .
cp ../04-NPT2/NPT2ed.gro .
cp ../00-Prep/posre*.itp .
cp ../00-Prep/index_$FN.ndx ./i.ndx
cp ../01-Min/lig*.itp .
cp ../01-Min/min.tpr .

$GMX  grompp -f md.mdp -c NPT2ed.gro -p $FN.top -o md.tpr -t NPT2.cpt -maxwarn 1 -r NPT2ed.gro -n i.ndx

mpirun -np ${SLURM_NTASKS_PER_NODE} --bind-to socket \ 
	mdrun_mpi -s md.tpr -deffnm MD -maxh 71.0  \
        -ntomp ${SLURM_CPUS_PER_TASK} &> run.out

echo 0 | $GMX trjconv -s md.tpr -f md.trr -o md_lastframe.gro -b 30000 -e 30000 -pbc whole
echo 0 | $GMX trjconv -s md.tpr -f md.trr -o md_reimaged.xtc -pbc whole

echo Protein Protein | $GMX trjconv -s min.tpr -f md_reimaged.xtc -o md_protraj.xtc -pbc nojump -center
echo Protein Protein | $GMX trjconv -s min.tpr -f md_protraj.xtc -o md_cluster.xtc -pbc cluster
echo Protein Protein | $GMX trjconv -s min.tpr -f md_cluster.xtc -o md_final.xtc -pbc mol -ur compact -center
echo Protein | $GMX gyrate -s min.tpr -f md_final.xtc -o ${FN}_Rgyr.xvg
echo Backbone Backbone | $GMX rms -s min.tpr -f md_final.xtc -o ${FN}_RMSD.xvg -fit rot+trans

echo Protein | $GMX trjconv -s min.tpr -f md_final.xtc -o ${FN}_protein.pdb -b 30000 -e 30000
echo Pressure | $GMX energy -f md.edr -o ${FN}_pressure.xvg

echo -e "\nStarting creation of pressure.dat file \n"
grep -v '[@#]' ${FN}_pressure.xvg | grep "  0.9" >> ${FN}_pressure1bar.dat
grep "  1.00" ${FN}_pressure.xvg >> ${FN}_pressure1bar.dat
grep "  1.01" ${FN}_pressure.xvg >> ${FN}_pressure1bar.dat
grep "  1.02" ${FN}_pressure.xvg >> ${FN}_pressure1bar.dat
echo -e "SUCCESS: pressure.dat file created\n"


