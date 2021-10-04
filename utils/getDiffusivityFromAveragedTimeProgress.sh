#!/bin/bash

fname=$1
ZE=$2


D=`awk -vZE=$ZE <"$fname" 'BEGIN {count=0; summ=0.0} NR>2 {if($2>ZE){summ+=(log($2**2/$1)); count+=1 } } END {print exp(summ/count)/2}'`

Derr=`awk -vZE=$ZE -vD=$D <"$fname" 'BEGIN {count=0; summ=0.0; ymean=log(2*D); } NR>2 {if($2>ZE){summ+=(log($2**2/$1)-ymean)**2; count+=1 } } END {ysem=sqrt(summ/count); print 0.5*ysem*D }'`


echo $D" "$Derr
