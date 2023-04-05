#!/bin/bash

name=FILE_STEM_HERE

sander -O -i min0.in -c ${name}.rst7 -p ${name}.prmtop -o ${name}_min0.out -r ${name}_min0.rst7

sander -O -i min1.in -c ${name}_min0.rst7 -p ${name}.prmtop -o ${name}_min1.out -r ${name}_min1.rst7

sander -O -i min2.in -c ${name}_min1.rst7 -p ${name}.prmtop -o ${name}_min2.out -r ${name}_min2.rst7
