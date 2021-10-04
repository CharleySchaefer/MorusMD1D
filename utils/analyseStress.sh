#!/bin/bash

script=$(readlink -f $0)
scriptpath=`dirname $script`
utils="$scriptpath"


if [ $# -lt 1 ] ; then
  echo "Usage: $0 <directory>"
  echo " "
  echo "This script reads param.in and timeprogress.out"
  echo "to calculate the stretch distribution."
  exit
fi
dir=$1
dir=${dir%/}
dir=${dir#./}
NBINS=100    # TODO: make input argument


# READ NUMBER OF ENTANGLEMENTS (ZE) FROM PARAM.IN
f_param=$dir/param.in
if test ! -f "$f_param"; then
    echo "error: $f_param does not exist."
    exit
fi
ZE=`awk < "$f_param" 'NR==1 {for (i=1; i<=NF; i++) { if ($(i)=="--Ze"){ print $(i+1); break}}} '`


# CHECK EXISTENCE OF TIMEPROGRESS.OUT
f_timeprogress=$dir/timeprogress.out
if test ! -f "$f_timeprogress"; then
    echo "error: $f_timeprogress does not exist."
    exit
fi


# GET STRESS DISTRIBUTION
mscript="addpath('$scriptpath');addpath('$utils');getStress('$f_timeprogress',$ZE,$NBINS)"
$utils/run_mscript.sh "$mscript" > $dir/StressDistribution.out
echo "Data written to $dir/StressDistribution.out"



