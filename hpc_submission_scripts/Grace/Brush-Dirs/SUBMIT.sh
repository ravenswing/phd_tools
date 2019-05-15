#!/bin/bash -l

cwd=$(pwd)
echo -e "\nCurrent working directory:"
echo -e "\t$cwd"

step=$(find . -maxdepth 1 -name '*.sh' -printf '%f\n' | grep -v SUBMIT.sh | grep -v run*.sh | grep -vi *cont.sh)
echo -e "\nScript file used to make run script:"
echo -e "\t$step"

echo "#!/bin/bash -l" > run.sh
echo "#$ -wd $cwd" >> run.sh
cat $step >> run.sh 

echo -e "\nSubmitting job to batch queue..."
sbatch run.sh
