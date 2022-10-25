#!/bin/bash

echo "#!/bin/bash" > _temp.sh

path=$(pwd)
echo -e "\nCurrent working directory:"
echo -e "\t$path"
echo "#SBATCH --workdir=$path" >> _temp.sh

step=$(find . -maxdepth 1 -name '*.sh' -printf '%f\n' | grep -v RUN.sh | grep -vi *T.sh | grep -vi _*.sh)
echo -e "\nScript file used to make run script:"
echo -e "\t$step"
cat $step >> _temp.sh

echo -e "\nSubmitting job to batch queue..."
sbatch _temp.sh
