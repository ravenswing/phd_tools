#!/bin/bash

# change the file name in the running script
export FN=$(basename -- "$PWD")

# replace run number in new script
sed -i "s/BASE_NAME/$FN/g" old_gpu_equil.sh

# submit job
qsub old_gpu_equil.sh

