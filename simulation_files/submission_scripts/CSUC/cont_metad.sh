 #!/bin/bash
#SBATCH -J PROD
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -N 1
#SBATCH -t 4-0
#SBATCH -e %x.%j.err
#SBATCH -o %x.%j.out


# --------------------------------- INPUTS ------------------------------------

steps_per_day=3		# Expected number of steps to run within walltime (24 hrs).
steps=50		# Maximum number of steps to run. 
traj=metad		# Name of all output files. 
ndx=i.ndx		# Name of index file.


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
#method=$(cd ../..; basename -- "$PWD")
plumed="plumed_${FN}.dat"

# Define the number steps and the max simulation time (ps).
tmax=$(($steps*10000))

# Put standard run information in the start of the log file.
printf "INFO | NODENAME:  %s \n" $SLURMD_NODENAME >> run.log
printf "INFO | RUNNING ON:  %s \n" $SLURM_SUBMIT_HOST >> run.log

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

for i in $(seq $start $stop)
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

# Relaunch if the total steps have not been reached.
if (( i < steps ))
then
   sbatch cont_metad.sh
else
    echo $(log_line "DONE" "All Steps ") >> run.log 
fi

exit
