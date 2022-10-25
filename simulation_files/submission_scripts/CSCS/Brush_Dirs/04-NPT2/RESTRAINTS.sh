#!/bin/bash

export FN=$(cd ..; basename -- "$PWD")

cp BACKUP.top $FN.top

head -n 14 posres_CAlpha.itp > C-alpha_pr.itp
    
sed -i -e '/EKEK_Protein/a \\n#ifdef CAPOSRES\n#include "C-alpha_pr.itp"\n#endif\n ' $FN.top

#for i in {2..12}
#do 
    #echo $i
    #sed -i -e '/EKEK_Protein${i}.itp/a \\n#ifdef CAPOSRES\n#include "C-alpha_pr.itp"\n#endif\n' $FN.top
#done


