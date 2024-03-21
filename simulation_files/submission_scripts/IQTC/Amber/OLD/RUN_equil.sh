#!/bin/bash

# change the file name in the running script
export FN=$(basename -- "$PWD")

# replace run number in new script
sed -i "s/BASE_NAME/$FN/g" iqtc_equil_old.sh

# submit job
qsub iqtc_equil_old.sh

