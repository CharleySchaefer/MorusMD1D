#!/bin/bash

# Sticker Friction: Experimental feature

# Simulation settings
Ze=200             # Number of entanglements per chain
SequenceN="0000000000"   # No stickers
SequenceS="0100000100"   # With stickers
Zs=2
Niter=5000000     # Number of time steps
Nprint=50000       # Print timeprogress after this number of steps
time_fac=0.4       # 
sticker_mobility=0.01

# Calculation of mu(t)
mu_N=50            # Calculate statistics after this number of steps
tL=1e-1            # mu(t) histogram: lower time boundary
tU=1e11            # mu(t) histogram: upper time boundary
nbins=100     # mu(t) histogram: number of time bins

fnameN="demo_sticker_frictionN.txt"
fnameS="demo_sticker_frictionS.txt"
#-----------------------------------------------------
# FIRST SIMULATION: WITH STICKERS
echo "Simulation with stickers"
./MorusMD --sequence $SequenceS --Ze $Ze --Ntime $Niter --Nprint $Nprint --mu-time-steps  $mu_N --mu-Nbins $nbins --mu-time-L $tL --mu-time-U $tU --time-fac $time_fac --sticker-mobility $sticker_mobility --sticker-all-open > $fnameS
# Analysis
Lmean=`awk <"timeprogress.out" 'BEGIN{count=0; summ=0.0} NR>1 {summ+=($3-$2); count++;} END {print summ/count}'`

Lstd=`awk -vLmean=$Lmean <"timeprogress.out" 'BEGIN{count=0; summ=0.0;} NR>1 {summ+=(($3-$2)-Lmean)**2; count++;} END {print sqrt(summ/count)}'`

Lsem=`awk -vLstd=$Lstd -vNsample=300 <"timeprogress.out" 'BEGIN {print Lstd/sqrt(Nsample-1)}'`

echo "Contour length"
echo " theoretical: 200+/-8.16"
echo " simulation:  "$Lmean"+/-"$Lstd" ; sem="$Lsem

#-----------------------------------------------------
# SECOND SIMULATION: WITHOUT STICKERS
echo " "
echo "Simulation without stickers"
./MorusMD --sequence $SequenceN --Ze $Ze --Ntime $Niter --Nprint $Nprint --mu-time-steps  $mu_N --mu-Nbins $nbins --mu-time-L $tL --mu-time-U $tU --time-fac $time_fac --sticker-mobility $sticker_mobility > $fnameN

# Analysis
Lmean=`awk <"timeprogress.out" 'BEGIN{count=0; summ=0.0} NR>1 {summ+=($3-$2); count++;} END {print summ/count}'`

Lstd=`awk -vLmean=$Lmean <"timeprogress.out" 'BEGIN{count=0; summ=0.0;} NR>1 {summ+=(($3-$2)-Lmean)**2; count++;} END {print sqrt(summ/count)}'`

Lsem=`awk -vLstd=$Lstd -vNsample=300 <"timeprogress.out" 'BEGIN {print Lstd/sqrt(Nsample-1)}'`

echo "Contour length"
echo " theoretical: 200+/-8.16"
echo " simulation:  "$Lmean"+/-"$Lstd" ; sem="$Lsem

#-----------------------------------------------------
# POSTPROCESS: Compare both simulations with theory
#              conclusion: improved agreement with theory
#                          for with increased Nbeads
../utils/Octave/get_mu_LMtheory.m $Ze $tL $tU $nbins > demo_LMtheory_Ze$Ze.txt


#-----------------------------------------------------
# PLOT

plot="set terminal pngcairo enhanced dashed
set output \"demo_sticker_friction.png\"

set xlabel \"t/(3Z_e^3)\"
set log x
set format x \"10^{%%T}\"

set ylabel \"{/Symbol m}(t)\"
set yrange [0:1]

set key left bottom 

set arrow from 1,0 to 1,1 nohead lw 1 dt 2 lc rgb \"#444444\" 

Zs=2
tauS=(1.0/$sticker_mobility)*($Ze/Zs)**2

plot \\
\"$fnameN\" u (\$1/(3*$Ze**3)):2           w p ps 1.3 pt 6 lc rgb \"#888888\" title \"N=10; Ze=$Ze, Zs=0\",\\
\"$fnameS\" u (\$1/(3*$Ze**3)):2           w p ps 1.3 pt 8 lc rgb \"#444444\" title \"N=10; Ze=$Ze, Zs=2\",\\
\"demo_LMtheory_Ze$Ze.txt\" u (\$1/(3*$Ze**3)):2 w l lw 2 lc rgb \"red\" title \"Likhtman-McLeish\" 
"

printf "$plot" | gnuplot
