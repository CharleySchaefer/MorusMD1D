#!/bin/bash

# Demo: 
#    > chain extension in flow without fluctuations
#      for 3 different boundary condition implementations (A, B and C ).
#    > chain extension with thermal fluctuations (boundary B) for various 
#       Weissenberg numbers + comparison to single-mode DEMG

# sequence: 11 nodes
# chain Ze=10 entanglements
# 3000 time steps
# export every 3 time steps
# flow rate: Weissenberg Wi=0.7
# export sequence

ZE=10 # fixed in all calculations

#==========================================================================
# PART i - No fluctuations ; vary boundary condition at same Weissenberg number
Wi=0.7
outdir=demo
mkdir -p $outdir
./MorusMD --outdir $outdir/A  --sequence 00000000000 --Ze $ZE \
           --set-boundary A --Ntime 3000 --Nprint 3 \
           --Wi $Wi --export-sequence --exclude-fluctuations
./MorusMD --outdir $outdir/B  --sequence 00000000000 --Ze $ZE \
           --set-boundary B --Ntime 3000 --Nprint 3 \
           --Wi $Wi --export-sequence --exclude-fluctuations
./MorusMD --outdir $outdir/C  --sequence 00000000000 --Ze $ZE \
           --set-boundary C --Ntime 3000 --Nprint 3 \
           --Wi $Wi --export-sequence --exclude-fluctuations

./MorusMD --outdir $outdir/A2  --sequence 00000000000 --Ze $ZE \
           --set-boundary A --Ntime 3000 --Nprint 3 \
           --Wi $Wi --export-sequence 
./MorusMD --outdir $outdir/B2  --sequence 00000000000 --Ze $ZE \
           --set-boundary B --Ntime 3000 --Nprint 3 \
           --Wi $Wi --export-sequence 
./MorusMD --outdir $outdir/C2  --sequence 00000000000 --Ze $ZE \
           --set-boundary C --Ntime 3000 --Nprint 3 \
           --Wi $Wi --export-sequence 


../utils/analyseStretch.sh $outdir/A
../utils/analyseStretch.sh $outdir/B
../utils/analyseStretch.sh $outdir/C
../utils/analyseStretch.sh $outdir/A2
../utils/analyseStretch.sh $outdir/B2
../utils/analyseStretch.sh $outdir/C2
../utils/analyseStress.sh $outdir/A
../utils/analyseStress.sh $outdir/B
../utils/analyseStress.sh $outdir/C
../utils/analyseStress.sh $outdir/A2
../utils/analyseStress.sh $outdir/B2
../utils/analyseStress.sh $outdir/C2


./MorusMD --outdir $outdir/L  --sequence 00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 --Ze $ZE \
           --set-boundary B --Ntime 2000000 --Nprint 2000 \
           --Wi $Wi  --exclude-fluctuations

for bnd in "A" "B" "C" ; do 
# Plot results
plot="set terminal pngcairo enhanced dashed
set output \"$outdir/demo_extensional_flow_Wi"$Wi"_boundary$bnd.png\"

set log x
set xlabel \"time/{/Symbol t}_e\"

set ylabel \"R/a\"

We=0.7
Z=10.0
pi=3.1415

Zstar=sin( pi*sqrt((We/Z**2))*(  Z/2  )  )/(pi*sqrt((We/Z**2))*cos( pi*sqrt((We/Z**2))*(Z/2)))

set style line 1 lt 1 lw 3 lc rgb \"red\"
set style line 2 lt 1 lw 1 lc rgb \"black\"
set key left bottom reverse

plot \\
\"$outdir/$bnd/timeprogress.out\" u 1:16  w l ls 2  title \"beads\", \\
\"\" u 1:6  w l ls 2 notitle, \\
\"\" u 1:7  w l ls 2 notitle, \\
\"\" u 1:8  w l ls 2 notitle, \\
\"\" u 1:9  w l ls 2 notitle, \\
\"\" u 1:10 w l ls 2 notitle, \\
\"\" u 1:11 w l ls 2 notitle, \\
\"\" u 1:12 w l ls 2 notitle, \\
\"\" u 1:13 w l ls 2 notitle, \\
\"\" u 1:14 w l ls 2 notitle, \\
\"\" u 1:15 w l ls 2 notitle, \\
\"\" u 1:2 w l ls 1 title \"end point\", \\
\"\" u 1:3  w l ls 1 notitle, \\
 Zstar w l lw 2 lt 2 lc rgb \"blue\" title \"analytical steady state\", \\
-Zstar w l lw 2 lt 2 lc rgb \"blue\" notitle
"
printf "$plot" | gnuplot
done # loop over boundart conditions finished
#=====================================================================================

#=====================================================================================
# Loop Weissenberg number with fluctuations

plotall="set terminal pngcairo enhanced
set output \"$outdir/demo_extensional_flow_demg.png\"

set log xy
set style line 1 lt 1 lw 2 lc rgb \"black\"
set style line 2 lt 1 lw 1 ps 0.8 pt 7 lc rgb \"#DD4444\"
set grid lt 1 lw 0.5 lc rgb \"grey\"

set xlabel \"t/{/Symbol t}_R\"
set ylabel \"{/Symbol l}\"
set yrange [0.5:1000]

plot 0 lw 0 notitle"
for Wi in  0.1 0.2 0.5 0.9 1.1 2 5 ; do
subdir=Wi$Wi
mkdir -p $outdir/$subdir
./MorusMD --outdir $outdir/$subdir --sequence 00000000000 --Ze $ZE --Ntime 3000 --Nprint 3 --Wi $Wi --export-sequence --time-fac 0.4

# Analysis
Lmean=`awk <"$outdir/$subdir/timeprogress.out" 'BEGIN{count=0; summ=0.0} NR>1 {summ+=($3-$2); count++;} END {print summ/count}'`
Lstd=`awk -vLmean=$Lmean <"$outdir/$subdir/timeprogress.out" 'BEGIN{count=0; summ=0.0;} NR>1 {summ+=(($3-$2)-Lmean)**2; count++;} END {print sqrt(summ/count)}'`
Lsem=`awk -vLstd=$Lstd -vNsample=300 <"$outdir/$subdir/timeprogress.out" 'BEGIN {print Lstd/sqrt(Nsample-1)}'`

#echo "Contour length"
#echo " theoretical: 10+/-1.826"
#echo " simulation:  "$Lmean"+/-"$Lstd" ; sem="$Lsem



# Plot results
plot="set terminal pngcairo enhanced dashed
set output \"$outdir/$subdir/demo_extensional_flow_Wi"$Wi"_2.png\"

set log x
set xlabel \"time/{/Symbol t}_e\"

set ylabel \"R/a\"

We=0.7
Z=10.0
pi=3.1415

Zstar=sin( pi*sqrt((We/Z**2))*(  Z/2  )  )/(pi*sqrt((We/Z**2))*cos( pi*sqrt((We/Z**2))*(Z/2)))

set style line 1 lt 1 lw 3 lc rgb \"red\"
set style line 2 lt 1 lw 1 lc rgb \"black\"
set key left bottom reverse

plot \\
\"$outdir/$subdir/timeprogress.out\" u 1:(\$16)  w l ls 2 title \"beads\", \\
\"\" u 1:(\$6)  w l ls 2 notitle, \\
\"\" u 1:(\$7)  w l ls 2 notitle, \\
\"\" u 1:(\$8)  w l ls 2 notitle, \\
\"\" u 1:(\$9)  w l ls 2 notitle, \\
\"\" u 1:(\$10) w l ls 2 notitle, \\
\"\" u 1:(\$11) w l ls 2 notitle, \\
\"\" u 1:(\$12) w l ls 2 notitle, \\
\"\" u 1:(\$13) w l ls 2 notitle, \\
\"\" u 1:(\$14) w l ls 2 notitle, \\
\"\" u 1:(\$15) w l ls 2 notitle, \\
\"\" u 1:(\$2) w l ls 1 title \"end point\", \\
\"\" u 1:(\$3)  w l ls 1 notitle, \\
 Zstar w l lw 2 lt 2 lc rgb \"blue\" title \"analytical steady state\", \\
-Zstar w l lw 2 lt 2 lc rgb \"blue\" notitle
"

printf "$plot" | gnuplot
# Plot results
plot="set terminal pngcairo enhanced dashed
set output \"$outdir/$subdir/demo_extensional_flow_Wi"$Wi"_3.png\"

set log x
set xlabel \"time/{/Symbol t}_e\"

set ylabel \"R/a\"

We=0.7
Z=10.0
pi=3.1415

Zstar=sin( pi*sqrt((We/Z**2))*(  Z/2  )  )/(pi*sqrt((We/Z**2))*cos( pi*sqrt((We/Z**2))*(Z/2)))

set style line 1 lt 1 lw 3 lc rgb \"red\"
set style line 2 lt 1 lw 1 lc rgb \"black\"
set key left bottom reverse

plot \\
\"$outdir/$subdir/timeprogress.out\" u 1:(\$16+\$4)  w l ls 2  title \"beads\", \\
\"\" u 1:(\$6+\$4)  w l ls 2 notitle, \\
\"\" u 1:(\$7+\$4)  w l ls 2 notitle, \\
\"\" u 1:(\$8+\$4)  w l ls 2 notitle, \\
\"\" u 1:(\$9+\$4)  w l ls 2 notitle, \\
\"\" u 1:(\$10+\$4) w l ls 2 notitle, \\
\"\" u 1:(\$11+\$4) w l ls 2 notitle, \\
\"\" u 1:(\$12+\$4) w l ls 2 notitle, \\
\"\" u 1:(\$13+\$4) w l ls 2 notitle, \\
\"\" u 1:(\$14+\$4) w l ls 2 notitle, \\
\"\" u 1:(\$15+\$4) w l ls 2 notitle, \\
\"\" u 1:(\$2+\$4)  w l ls 1 title \"end point\", \\
\"\" u 1:(\$3+\$4)  w l ls 1 notitle, \\
 Zstar w l lw 2 lt 2 lc rgb \"blue\" title \"analytical steady state\", \\
-Zstar w l lw 2 lt 2 lc rgb \"blue\" notitle
"

printf "$plot" | gnuplot
plotall="$plotall ,\
	\"demo/Wi$Wi/timeprogress.out\" u (\$1/($ZE*$ZE)):( (\$3-\$2)/$ZE ) w p ls 2 notitle,\
(exp( ($Wi-1)*x  )*$Wi-1)/($Wi-1) w l ls 1 notitle"
done # end loop weissenberg numbers


printf "$plotall" | gnuplot 
