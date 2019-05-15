#!/bin/bash -l

cwd=$(pwd)
echo -e "\nCurrent working directory:"
echo -e "\t$cwd"

export FN=$(cd ..; basename -- "$PWD")
num=$((`ls -ltrh *.cpt | grep -v prev | grep '_' | tail -n 1 | awk '{split($NF, a, "[_.]"); print a[length(a)-1]}'`+1))
echo -e "\nSystem and job number:"
echo -e "\t$FN\t$num"

echo "#!/bin/bash -l" > run2.sh
cat mdcont.sh >> run2.sh 
sed -i -e "/-pe mpi/a #$ -wd $cwd\n" run2.sh
sed -i -e "/-pe mpi/a #$ -N $FN-$num" run2.sh
sed -i "s/RUN_NUMBER/$num/g" run2.sh

echo -e "\nSubmitting job to batch queue..."
qsub run2.sh

