#!/bin/bash
#SBATCH -J PROD
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -N 1
#SBATCH -t 4-0
#SBATCH -e %x.%j.err
#SBATCH -o %x.%j.out


# --------------------------------- INPUTS ------------------------------------

steps=50		# Maximum number of steps to run. 
traj=metad		# Name of all output files. 
ndx=i.ndx		# Name of index file.
prev=md			# Name of previous .gro and .cpt files.

# --------------------------------- SETUP  ------------------------------------

# Load the modules and set the environment variables.
module load apps/gromacs/2023_plumed_2.9
export GMX=gmx_mpi
export GMX_DISABLE_GPU_TIMING=yes

# Create a custom logging function.
log_line() {
	local timestamp="$(date +"%F %T")"
	local line='.............................................'
	echo $(printf "%s  %s%s %s\n" "$timestamp" "$2" "${line:${#2}}" "$1")
}

# Get file names from the directory tree.
FN=$(cd ..; basename -- "$PWD")
# method=$(cd ../..; basename -- "$PWD")
plumed=plumed_${FN}.dat

# Define the number steps and the max simulation time (ps).
tmax=$(($steps*10000))

# Put standard run information in the start of the log file.
printf "INFO | SUBMIT_DIR:  %s \n" $SLURM_SUBMIT_DIR > run.log
printf "INFO | NODENAME:  %s \n" $SLURMD_NODENAME >> run.log
printf "INFO | RUNNING ON:  %s \n" $SLURM_SUBMIT_HOST >> run.log


# ----------------------------- INITIAL STEP  ---------------------------------

i=1

# Run grompp to initialise the simulation
$GMX grompp -f step10.mdp -c $prev.gro -p $FN.top -o ${traj}${i}.tpr -t $prev.cpt -r $prev.gro -n $ndx
echo $(log_line "RUNNING" "Step $i ") >> run.log

# Run Step 1.
OMP_NUM_THREADS=24 srun -n 1 -c 24 $GMX mdrun -dlb auto -pin auto       \
    -s ${traj}${i}.tpr -deffnm ${traj} -plumed $plumed

echo $(log_line "COMPLETED" "Step $i ") >> run.log

# Backup the checkpoint file and update the log file.
cp $traj.cpt ${traj}_step${i}.cpt


# ----------------------------- RUN NEXT STEPS --------------------------------

for i in $(seq 2 $steps)
do	
    # Extend the previous tpr by 10 ns.
    $GMX convert-tpr -s ${traj}$((i-1)).tpr -o ${traj}${i}.tpr -extend 10000
    echo $(log_line "RUNNING" "Step $i ") >> run.log
 
    # Run the next step. 
    OMP_NUM_THREADS=24 srun -n 1 -c 24 $GMX mdrun -dlb auto -pin auto   \
        -s ${traj}${i}.tpr -deffnm ${traj} -plumed $plumed              \
        -cpi ${traj}_step$((i-1)).cpt -noappend 
    
    echo $(log_line "COMPLETED" "Step $i ") >> run.log
    
    # Backup checkpoint file and update the log file.
    cp $traj.cpt ${traj}_step${i}.cpt 
done


# ------------------------------ FINALISING  ----------------------------------

exit
