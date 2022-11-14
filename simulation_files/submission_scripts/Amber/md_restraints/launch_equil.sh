#!/bin/bash

name="ship1"

pmemd.cuda -O -i eq_1.in -p ${name}.top -c ${name}_min2.rst7 -r ${name}.eq_1.r -x ${name}.eq_1.x -e ${name}.eq_1.e -o ${name}.eq_1.o -ref ${name}_min2.rst7 -inf eq_1.inf
gzip -f ${name}.eq_1.e ${name}.eq_1.o

pmemd.cuda -O -i eq_2.in -p ${name}.top -c ${name}.eq_1.r -r ${name}.eq_2.r -x ${name}.eq_2.x -e ${name}.eq_2.e -o ${name}.eq_2.o -ref ${name}.eq_1.r -inf eq_2.inf
gzip -f ${name}.eq_2.e ${name}.eq_2.o

pmemd.cuda -O -i eq_3.in -p ${name}.top -c ${name}.eq_2.r -r ${name}.eq_3.r -x ${name}.eq_3.x -e ${name}.eq_3.e -o ${name}.eq_3.o -ref ${name}.eq_2.r -inf eq_3.inf
gzip -f ${name}.eq_3.e ${name}.eq_3.o

pmemd.cuda -O -i eq_4.in -p ${name}.top -c ${name}.eq_3.r -r ${name}.eq_4.r -x ${name}.eq_4.x -e ${name}.eq_4.e -o ${name}.eq_4.o -ref ${name}.eq_3.r -inf eq_4.inf
gzip -f ${name}.eq_4.e ${name}.eq_4.o

pmemd.cuda -O -i eq_5.in -p ${name}.top -c ${name}.eq_4.r -r ${name}.eq_5.r -x ${name}.eq_5.x -e ${name}.eq_5.e -o ${name}.eq_5.o -ref ${name}.eq_4.r -inf eq_5.inf
gzip -f ${name}.eq_5.e ${name}.eq_5.o

pmemd.cuda -O -i eq_6.in -p ${name}.top -c ${name}.eq_5.r -r ${name}.eq_6.r -x ${name}.eq_6.x -e ${name}.eq_6.e -o ${name}.eq_6.o -ref ${name}_eq_5.r -inf eq_6.inf
gzip -f ${name}.eq_6.e ${name}.eq_6.o

exit
