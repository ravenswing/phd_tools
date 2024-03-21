#!/bin/csh 
#$ -S /bin/bash
#$ -N Back-Funnel-a2b1-SC4
#$ -pe smp 1
#$ -q iqtc10_g1819_gpu.q
#$ -cwd 
#$ -o $JOB_NAME.out
#$ -e $JOB_NAME.err

# --------------------------------- INPUTS ------------------------------------

steps_per_day=3		# Expected number of steps to run within walltime (24 hrs).
steps=50		# Maximum number of steps to run. 
traj=metad		# Name of all output files. 
ndx=i.ndx		# Name of index file.


# --------------------------------- SETUP  ------------------------------------

# Load the modules and set the environment variables.
source /etc/profile
module load gromacs/2021.4_cuda
export OMP_NUM_THREADS=1
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
method=$(cd ../..; basename -- "$PWD")
plumed=${FN}_${method}.dat

# Define the number steps and the max simulation time (ps).
tmax=$(($steps*10000))

# Put standard run information in the start of the log file.
printf "INFO | TMPDIR:  %s \n" $TMPDIR >> run.log

# Get the node information for error tracking.
rsync -a run.log $TMPDIR/
cd $TMPDIR 
printf "INFO | RUNNING ON:  %s \n" $HOSTNAME >> run.log
rsync -a run.log $SGE_O_WORKDIR/.
cd $SGE_O_WORKDIR

# Calculate the previously completed step (based on name, NOT time created).
start=$((`ls -lv ${traj}*.cpt | grep '.' | tail -n 1 | awk '{split($NF, a, "[_.]"); print a[length(a)-1]}' | awk '{split($NF, b, "step"); print b[length(b)]}'`+1))
# Set the numbers of steps to run in the walltime...
stop=$(( $start + $steps_per_day - 1))
# ... and adjust if it will exceed max. required steps. 
if (( stop > steps ))
then
   stop=$steps
fi

echo $(log_line "RESTARTING" "From Step $start ==> $stop ") >> run.log

# ----------------------------- RUN NEXT STEPS --------------------------------

# Copy all metaD files -> node directory.
rsync -a $plumed $TMPDIR/
rsync -a *.pdb $TMPDIR/
rsync -a COLVAR $TMPDIR/
rsync -a HILLS $TMPDIR/

for i in $(seq $start $stop)
do	
    # Extend the previous tpr by 10 ns.
    $GMX convert-tpr -s ${traj}$((i-1)).tpr -o ${traj}${i}.tpr -extend 10000

    # Copy step inputs to the directory where the jobs will run.
    rsync -au ${traj}${i}.tpr $TMPDIR/
    rsync -au ${traj}_step$((i-1)).cpt $TMPDIR/
 
    echo $(log_line "RUNNING" "Step $i ") >> run.log
    rsync -a run.log $TMPDIR/
    
    # Run the next step.
    cd $TMPDIR 
    $GMX mdrun -s ${traj}${i}.tpr -deffnm ${traj} -cpi ${traj}_step$((i-1)).cpt -noappend -plumed $plumed
    echo $(log_line "COMPLETED" "Step $i ") >> run.log
    rsync -au * $SGE_O_WORKDIR/.
    cd $SGE_O_WORKDIR
    
    # Backup checkpoint file and update the log file.
    cp $traj.cpt ${traj}_step${i}.cpt 
    echo $(log_line "SYNCED" "Step $i ") >> run.log
done

# ------------------------------ FINALISING  ----------------------------------

# Account for system clock problems
touch *

# Relaunch if the total steps have not been reached.
if (( i < steps ))
then
   qsub cont_metad.sh
else
    echo $(log_line "DONE" "All Steps ") >> run.log 
fi

exit
