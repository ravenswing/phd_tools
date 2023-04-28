#!/bin/bash -l

# open the new run script
echo "#!/bin/bash -l" > _run_ext.sh

# add the current working dir.
path=$(pwd)
echo -e "\nCurrent working directory:"
echo -e "\t$path"
echo "#SBATCH --chdir=$path" >> _run_ext.sh

# calculate current run number
num=$((`ls -ltrh *.cpt | grep -v prev | grep '_' | tail -n 1 | awk '{split($NF, a, "[_.]"); print a[length(a)-1]}'`+1))

# edit job-name
export FN=$(cd ..; basename -- "$PWD")
echo "#SBATCH --job-name=$FN-$num" >> _run_ext.sh

# import rest of the commands from cont file
cat ext.sh >> _run_ext.sh 

# replace run number in new script
sed -i "s/RUN_NUMBER/$num/g" _run_ext.sh

# submit job
#sbatch _run_ext.sh
