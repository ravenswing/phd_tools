#!/bin/csh 
#$ -S /bin/bash
#$ -N POSTPROC
#$ -pe smp 1
#$ -q iqtc10_g1819_gpu.q
#$ -cwd 
#$ -o $JOB_NAME.out
#$ -e $JOB_NAME.err

# --------------------------------- INPUTS ------------------------------------

steps=10		# Maximum number of steps run.
traj=md			# Name of all output files. 
ndx=i.ndx		# Name of index file.
tpr=md9.tpr		# Tpr to use for trjconv transformations.

# --------------------------------- SETUP  ------------------------------------

# Load the modules and set the environment variables.
source /etc/profile
module load gromacs/2021.4_cuda
export OMP_NUM_THREADS=24
export SGE_ROOT="/sge"
export GMX=gmx_mpi
export CUDA_VISIBLE_DEVICES=`cat $TMPDIR/.gpus`

# Create a custom logging function.
log_line() {
	local timestamp="$(date +"%F %T")"
	local line='.............................................'
	echo $(printf "%s  %s%s %s\n" "$timestamp" "$2" "${line:${#2}}" "$1")
}

# Get file names from the directory tree.
FN=$(cd ..; basename -- "$PWD")

# Define the number steps and the max simulation time (ps).
tmax=$(($steps*100))

# Put standard run information in the start of the log file.
printf "INFO | TMPDIR:  %s \n" $TMPDIR >> run.log

# Get the node information for error tracking.
rsync -a run.log $TMPDIR/
cd $TMPDIR 
printf "INFO | RUNNING ON:  %s \n" $HOSTNAME >> run.log
rsync -a run.log $SGE_O_WORKDIR/.
cd $SGE_O_WORKDIR

# ------------------------------- POST-PROCESSING ----------------------------------

out_name=$FN

echo $(log_line "RUNNING" "Post-processing ") >> run.log

# Concatenate the trajectories.
$GMX trjcat -f ./$traj*.trr -o $out_name.trr
$GMX trjcat -f ./$traj*.xtc -o $out_name.xtc

# Reformat the trajectories. 
echo Protein_LIG         | $GMX trjconv -s $tpr -f ${out_name}.xtc       -o ${out_name}_whole.xtc -pbc whole -n i.ndx
echo Protein Protein_LIG | $GMX trjconv -s $tpr -f ${out_name}_whole.xtc -o ${out_name}_clust.xtc -pbc cluster -n i.ndx
echo Protein Protein_LIG | $GMX trjconv -s $tpr -f ${out_name}_clust.xtc -o ${out_name}_final.xtc -pbc mol -ur compact -center -n i.ndx
# Save the final frame as a dry pdb. 
echo Protein_LIG         | $GMX trjconv -s $tpr -f ${out_name}_final.xtc -o ${out_name}_lastframe.pdb -b $tmax -e $tmax -n i.ndx

echo $(log_line "COMPLETED" "Post-processing ") >> run.log

exit
