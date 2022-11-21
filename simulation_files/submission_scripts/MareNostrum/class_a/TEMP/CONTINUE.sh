#!/bin/bash -l

# open the new run script
echo "#!/bin/bash -l" > _run3.sh

# add the current working dir.
path=$(pwd)
echo -e "\nCurrent working directory:"
echo -e "\t$path"
echo "#SBATCH --chdir=$path" >> _run3.sh

# calculate current run number
num=$((`ls -ltrh *.cpt | grep -v prev | grep '_' | tail -n 1 | awk '{split($NF, a, "[_.]"); print a[length(a)-1]}'`+1))

# edit job-name
export FN=$(cd ..; basename -- "$PWD")
echo "#SBATCH --job-name=$FN-$num" >> _run3.sh

# import rest of the commands from cont file
cat cont.sh >> _run3.sh 

# replace run number in new script
sed -i "s/RUN_NUMBER/$num/g" _run3.sh

# submit job
sbatch _run3.sh
