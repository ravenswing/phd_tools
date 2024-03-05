#!/bin/bash

export name=$(cd ..; basename -- "$PWD")
export GMX=gmx_mpi

# copy relevant gmx files from prep dir.
#cp ../00-Prep/$name.top ./
#cp ../00-Prep/$name.gro ./
#cp ../00-Prep/*.itp ./
#cp ../00-Prep/i.ndx .

# run the minimisation
$GMX grompp -f min.mdp -c $name.gro -p $name.top -o min.tpr -n i.ndx
$GMX mdrun -s min.tpr -deffnm min -ntomp 8 -pin on

#Monitor the energy
echo 10 0 | $GMX energy -f min.edr -o ${name}_min_energy.xvg
#produce the gro file for NVT equilibration
echo 0 | $GMX trjconv -s min.tpr -f min.trr -o min_out.gro -pbc whole
#produce a readable minimised pdb for reference
echo 0 | $GMX trjconv -s min.tpr -f min.trr -o ${name}_minimised.pdb -pbc mol -ur compact
echo Protein | $GMX trjconv -s min.tpr -f min.trr -o ${name}_minimised_dry.pdb -pbc mol -ur compact -n i.ndx

