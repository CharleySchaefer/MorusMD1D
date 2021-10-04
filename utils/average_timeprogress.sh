#!/bin/bash


script=$(readlink -f $0)
scriptpath=`dirname $script`
utils=$scriptpath

dir=$1

$utils/run_mscript.sh "addpath('$utils');average_timeprogress('$dir')"


