#!/bin/bash -l

echo "#!/bin/bash -l" > _run.sh

path=$(pwd)
echo -e "\nCurrent working directory:"
echo -e "\t$path"
echo "#SBATCH --chdir=$path" >> _run.sh

name=$(cd ..; basename -- "$PWD")
echo "#SBATCH --job-name=$name0" >> _run.sh

cat prod.sh >> _run.sh 

echo -e "\nSubmitting job to batch queue..."
sbatch _run.sh
