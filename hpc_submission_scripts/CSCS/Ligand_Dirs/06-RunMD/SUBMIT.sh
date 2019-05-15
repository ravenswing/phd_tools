#!/bin/bash -l

export FN=$(cd ..; basename -- "$PWD")
echo "#!/bin/bash -l" > run.sh
echo "#SBATCH --job-name=MD_$FN" >> run.sh
cat mdrun.sh >> run.sh 

sbatch run.sh
