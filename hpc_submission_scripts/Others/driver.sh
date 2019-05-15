#!/bin/bash
#SBATCH --job-name="DRVR_MOLECULE" 
#SBATCH --time=25:00:00
#SBATCH --error=MDrun_errors.%j
#SBATCH --output=MDrun_output.%j
#SBATCH --nodes=1

module load plumed/2.3.3_libmatheval

export FN=$(basename -- "$PWD")

#plumed driver --plumed plumed.dat --mf_xtc MDrun_final.xtc --timestep 0.002
plumed driver --plumed plumed.dat --mf_xtc md_final.xtc --timestep 0.002

