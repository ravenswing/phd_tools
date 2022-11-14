#!/bin/bash

name="FILE_STEM_HERE"

echo "Starting Minimisation 0"
sander -O -i min0.in -c ${name}.crd -p ${name}.top -o ${name}_min0.out -r ${name}_min0.rst

echo "Starting Minimisation 1"
sander -O -i min1.in -c ${name}_min0.rst -p ${name}.top -o ${name}_min1.out -r ${name}_min1.rst

echo "Starting Minimisation 2"
sander -O -i min2.in -c ${name}_min1.rst -p ${name}.top -o ${name}_min2.out -r ${name}_min2.rst
