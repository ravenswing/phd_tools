#!/bin/bash -l

echo "#!/bin/bash" > _run.sh

#path=$(pwd)
#echo -e "\nCurrent working directory:"
#echo -e "\t$path"
#echo "#SBATCH --chdir=$path" >> _run.sh

name1=$(cd ../..; basename -- "$PWD")
name2=$(cd ..; basename -- "$PWD")
name3=$(basename -- "$PWD")
echo "#SBATCH -J $name1-$name2-$name3" >> _run.sh

cat benchmark.sh >> _run.sh 

echo -e "\nSubmitting job to batch queue..."
sbatch _run.sh
