#!/bin/csh 
#$ -S /bin/bash
#$ -N Back-Funnel-a2b1-SC4
#$ -pe smp 1
#$ -q iqtc10_g1819_gpu.q
#$ -cwd 
#$ -o $JOB_NAME.out
#$ -e $JOB_NAME.err

# --------------------------------- INPUTS ------------------------------------

steps=50		# Maximum number of steps to run. 
traj=metad		# Name of all output files. 
ndx=i.ndx		# Name of index file.
prev=md			# Name of previous .gro and .cpt files.

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
printf "INFO | SGE_O_WORKDIR:  %s \n" $SGE_O_WORKDIR > run.log
printf "INFO | NSLOTS:  %s \n" $NSLOTS >> run.log
printf "INFO | TMPDIR:  %s \n" $TMPDIR >> run.log

# Get the node information for error tracking.
rsync -a run.log $TMPDIR/
cd $TMPDIR 
printf "INFO | RUNNING ON:  %s \n" $HOSTNAME >> run.log
rsync -a run.log $SGE_O_WORKDIR/.
cd $SGE_O_WORKDIR


# ----------------------------- INITIAL STEP  ---------------------------------

i=1

# Run grompp to initialise the simulation
$GMX grompp -f step10.mdp -c $prev.gro -p $FN.top -o ${traj}${i}.tpr -t $prev.cpt -r $prev.gro -n $ndx

# Copy inputs and files needed to the directory where the jobs will run.
rsync -au ${traj}${i}.tpr $TMPDIR/
rsync -au $plumed $TMPDIR/
rsync -au *.pdb $TMPDIR/

echo $(log_line "RUNNING" "Step $i ") >> run.log
rsync -a run.log $TMPDIR/

# Run Step 1 on the compute node.
cd $TMPDIR 
$GMX mdrun -s ${traj}${i}.tpr -deffnm ${traj} -plumed $plumed
echo $(log_line "COMPLETED" "Step $i ") >> run.log
rsync -au * $SGE_O_WORKDIR/.
cd $SGE_O_WORKDIR

# Backup the checkpoint file and update the log file.
cp $traj.cpt ${traj}_step${i}.cpt
echo $(log_line "SYNCED" "Step $i ") >> run.log

# ----------------------------- RUN NEXT STEPS --------------------------------

for i in $(seq 2 $steps)
do	
    # Extend the previous tpr by 10 ns.
    $GMX convert-tpr -s ${traj}$((i-1)).tpr -o ${traj}${i}.tpr -extend 10000

    # Copy inputs and files needed to the directory where the jobs will run.
    rsync -a ${traj}${i}.tpr $TMPDIR/
    rsync -a ${traj}_step$((i-1)).cpt $TMPDIR/
    # rsync -a COLVAR $TMPDIR/
    # rsync -a HILLS $TMPDIR/
 
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

exit
