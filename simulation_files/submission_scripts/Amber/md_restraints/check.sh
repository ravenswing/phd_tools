#!/bin/bash

rm -f -- CHECK

for file in `ls -tr *md_*.x`;
do
	num=$(echo "$file" | awk -F'_' '{print $NF}' | awk -F'.' '{print $1}')
	frames=$(ncdump -h $file | grep 'frame =' | awk -F'(' '{print $NF}' | awk -F')' '{print $1}')

	echo "Step: $num    Frames: $frames" >> CHECK    
done

exit
