#!/bin/bash -l

echo "#!/bin/bash" > _run.sh

#path=$(pwd)
#echo -e "\nCurrent working directory:"
#echo -e "\t$path"
#echo "#SBATCH --chdir=$path" >> _run.sh

name=$(cd ..; basename -- "$PWD")
name2=$(basename -- "$PWD")
echo "#SBATCH -J $name-$name2" >> _run.sh

cat ext.sh >> _run.sh 

echo -e "\nSubmitting job to batch queue..."
sbatch _run.sh
