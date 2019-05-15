#!/bin/bash -l

export FN=$(cd ..; basename -- "$PWD")
num=$((`ls -ltrh *.cpt | grep -v prev | grep '_' | tail -n 1 | awk '{split($NF, a, "[_.]"); print a[length(a)-1]}'`+1))

echo "#!/bin/bash -l" > run2.sh
echo "#SBATCH --job-name=$FN-$num" >> run2.sh
cat mdcont.sh >> run2.sh 
sed -i "s/RUN_NUMBER/$num/g" run2.sh

sbatch run2.sh
