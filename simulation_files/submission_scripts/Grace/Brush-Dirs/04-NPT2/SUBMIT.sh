#!/bin/bash -l

cwd=$(pwd)
echo -e "\nCurrent working directory:"
echo -e "\t$cwd"

step=$(find . -maxdepth 1 -name '*.sh' -printf '%f\n' | grep -v SUBMIT.sh | grep -v run*.sh | grep -vi *cont.sh)
echo -e "\nScript file used to make run script:"
echo -e "\t$step"

echo "#!/bin/bash -l" > run.sh
cat $step >> run.sh 
sed -i -e "/-pe mpi/a #$ -wd $cwd\n" run.sh

echo -e "\nSubmitting job to batch queue..."
qsub run.sh
