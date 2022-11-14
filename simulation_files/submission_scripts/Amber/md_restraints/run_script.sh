#!/bin/bash

name="ship1"

pmemd.cuda -O -i mdin -p ${name}.top -c ${name}.eq_6.r -r ${name}.md_1.r -x ${name}.md_1.x -e ${name}.md_1.e -o ${name}.md_1.o -ref ${name}.eq_6.r -inf md_1.inf
gzip -f ${name}.md_1.e ${name}.md_1.o

for i in {2..200..1}
	do
		echo "Running Step: $i" >> run.log
		pmemd.cuda -O -i mdin -p ${name}.top -c ${name}.md_$((i-1)).r -r ${name}.md_$i.r -x ${name}.md_$i.x -e ${name}.md_$i.e -o ${name}.md_$i.o -ref  ${name}.md_$((i-1)).r -inf md_$i.inf
		gzip -f ${name}.md_$i.e ${name}.md_$i.o
	done
exit
