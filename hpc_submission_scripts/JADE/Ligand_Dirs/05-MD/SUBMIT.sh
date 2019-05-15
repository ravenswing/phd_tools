#!/bin/bash -l

export FN=$(cd ..; basename -- "$PWD")
echo "#!/bin/bash" > run.sh
echo "#SBATCH -J MD_$FN" >> run.sh
cat md.sh >> run.sh 

sbatch run.sh
