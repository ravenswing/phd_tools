#!/bin/bash -l

echo "#!/bin/bash" > _run.sh

#path=$(pwd)
#echo -e "\nCurrent working directory:"
#echo -e "\t$path"
#echo "#SBATCH --chdir=$path" >> _run.sh

# calculate current run number        
num=$((`ls -ltrh *.cpt | grep -v prev | grep '_' | tail -n 1 | awk '{split($NF, a, "[_.]"); print a[length(a)-1]}'`+1))

name=$(cd ..; basename -- "$PWD")
name2=$(basename -- "$PWD")
echo "#SBATCH -J $name-$name2" >> _run.sh

cat cont.sh >> _run.sh 

echo -e "\nSubmitting job to batch queue..."
sbatch _run.sh
