#!/bin/bash -l

# open the new run script
echo "#!/bin/bash -l" > _post.sub

# add the current working dir.
path=$(pwd)
echo -e "\nCurrent working directory:"
echo -e "\t$path"
echo "#SBATCH --chdir=$path" >> _post.sub

# edit job-name
name=$(cd ..; basename -- "$PWD")
name2=$(basename -- "$PWD")
echo "#SBATCH --job-name=pp$name-$name2" >> _post.sub

# import rest of the commands from cont file
cat postproc.sh >> _post.sub 

# replace run number in new script
sed -i "s/RUN_NUMBER/$num/g" _post.sub

# submit job
sbatch _post.sub
