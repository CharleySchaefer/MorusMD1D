#!/bin/bash

fname=$1 # timeprogress file
fname=${fname#./}
fname=${fname%/}
Nheader=2
script=$(readlink -f $0)
scriptpath=`dirname $script`


Lmean=`awk -vNheader=$Nheader <"$fname" 'BEGIN{count=0; summ=0.0} NR>Nheader {summ+=($3-$2); count++;} END {print summ/count}'`

Lstd=`awk -vLmean=$Lmean -vNheader=$Nheader<"$fname" 'BEGIN{count=0; summ=0.0;} NR>Nheader {summ+=(($3-$2)-Lmean)**2; count++;} END {print sqrt(summ/count)}'`

printf "%16s %12s\n" "#ContourLength" "Fluctuation"
printf "%16s %12s\n" $Lmean $Lstd
