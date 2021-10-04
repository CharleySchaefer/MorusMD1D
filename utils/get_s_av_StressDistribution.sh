#!/bin/bash

dir=$1

script=$(readlink -f $0)
scriptpath=`dirname $script`
utils="$scriptpath"

Nseeds=0
for seed in `find $dir/seed* -maxdepth 0` ; do
  Nseeds=$((Nseeds + 1))
done


mscript="
  function av_stress()
  Nseeds=$Nseeds;
  pth='$dir';
    for i=1:Nseeds
      fname=sprintf('%%s/seed%%d/StressDistribution.out', pth, i);
      data=importdata(fname);
      if i==1
        X=data(:,1);
        Y=data(:,2);
      else
        Y=data(:,2)+Y;
      end
    end
    Y=Y/Nseeds; 

    YVAR=zeros(size(Y));
    for i=1:Nseeds
      fname=sprintf('%%s/seed%%d/StressDistribution.out', pth, i);
      data=importdata(fname);
      YVAR=YVAR+(Y-data(:,2)).^2;
    end
    YVAR=YVAR/Nseeds; 

    for i=1:length(YVAR)
      fprintf('%%12e %%12e %%12e %s', X(i), Y(i), sqrt(YVAR(i))/(Nseeds-1)  );
    end
    end
" 
printf "$mscript" \\n > $dir/av_stress.m
$utils/run_mscript.sh "addpath('$dir');av_stress" > $dir/StressDistribution.out

