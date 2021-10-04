#!/bin/bash
Ze=10
tL=1e-2
tU=1e5
nbins=200

dir=demo/LMtheory
mkdir -p demo
mkdir -p $dir 
do_simN3=0 ; fnameN3="$dir/mu_Ze10_N3.out"
do_simN5=0 ; fnameN5="$dir/mu_Ze10_N5.out"
do_simN10=0 ; fnameN10="$dir/mu_Ze10_N10.out"
do_simN20=0 ; fnameN20="$dir/mu_Ze10_N20.out"
do_simN50=0 ; fnameN50="$dir/mu_Ze10_N50.out"
do_simN100=0 ; fnameN100="$dir/mu_Ze10_N100.out"
do_simN200=0 ; fnameN200="$dir/mu_Ze10_N200.out"
  
if [ $do_simN3 -eq 1 ] ; then
./MorusMD --set-boundary A --Ze $Ze --sequence 000 --export-mu --mu-time-steps 1 --mu-time-L $tL --mu-time-U $tU --Ntime 100000000 --time-fac 0.01 > $fnameN3
fi
  
if [ $do_simN5 -eq 1 ] ; then
./MorusMD --set-boundary A --Ze $Ze --sequence 00000 --export-mu --mu-time-steps 1 --mu-time-L $tL --mu-time-U $tU --Ntime 100000000 --time-fac 0.01 > $fnameN5
fi
  
if [ $do_simN10 -eq 1 ] ; then
./MorusMD --set-boundary A --Ze $Ze --sequence 0000000000 --export-mu --mu-time-steps 1 --mu-time-L $tL --mu-time-U $tU --Ntime 100000000 --time-fac 0.01 > $fnameN10
fi
  
if [ $do_simN20 -eq 1 ] ; then
./MorusMD --set-boundary A --Ze $Ze --sequence 00000000000000000000 --export-mu --mu-time-steps 1 --mu-time-L $tL --mu-time-U $tU --Ntime 100000000 --time-fac 0.05 > $fnameN20
fi
  
if [ $do_simN50 -eq 1 ] ; then
./MorusMD --set-boundary A --Ze $Ze --sequence 00000000000000000000000000000000000000000000000000 --export-mu --mu-time-steps 1 --mu-time-L $tL --mu-time-U $tU --Ntime 100000000 --time-fac 0.05 > $fnameN50
fi
  
if [ $do_simN100 -eq 1 ] ; then
./MorusMD --set-boundary A --Ze $Ze --sequence 0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 --export-mu --mu-time-steps 1 --mu-time-L $tL --mu-time-U $tU --Ntime 100000000 --time-fac 0.1 > $fnameN100
fi
  
if [ $do_simN200 -eq 1 ] ; then
./MorusMD --set-boundary A --Ze $Ze --sequence 00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 --export-mu --mu-time-steps 1 --mu-time-L $tL --mu-time-U $tU --Ntime 100000000 --time-fac 0.1 > $fnameN200
fi

#-----------------------------------------------------
# POSTPROCESS: Compare both simulations with theory
#              conclusion: improved agreement with theory
#                          for with increased Nbeads
utils="../utils"
$utils/run_mscript.sh "addpath('$utils/Octave') ; get_mu_LMtheory($Ze, $tL, $tU, $nbins); " > $dir/LMtheory_Ze$Ze.txt



plot="set terminal pngcairo enhanced
set output \"$dir/demo_LM_Ze10.png\"

set xlabel \"t/{/Symbol t}_e\"
set log x
set xtics format \"10^{%%T}\"

set ylabel \"{/Symbol m}(t)\"
set yrange [0:1]

set label \"Ze=10\" at  7000,0.9

set grid lw 0.5 lt 1 lc rgb \"grey\" 

set key left bottom reverse Left

plot \\
\"$fnameN3\"   u (\$1):2 w p ps 1.1 pt 4 lc rgb \"#CCCCCC\"  title \"N=3\",\
\"$fnameN5\"   u (\$1):2 w p ps 1.1 pt 6 lc rgb \"#AAAAAA\"  title \"N=5\",\
\"$fnameN10\"  u (\$1):2 w p ps 1.0 pt 8 lc rgb \"#888888\"  title \"N=10\",\
\"$fnameN20\"  u (\$1):2 w p ps 1.0 pt 10 lc rgb \"#666666\" title \"N=50\",\
\"$fnameN50\"  u (\$1):2 w p ps 0.9 pt 12 lc rgb \"#444444\" title \"N=50\",\
\"$fnameN100\" u (\$1):2 w p ps 0.9 pt 14 lc rgb \"#222222\" title \"N=100\",\
\"$fnameN200\" u (\$1):2 w p ps 0.5 pt 7 lc rgb \"black\" title \"N=200\",\
\"$dir/LMtheory_Ze$Ze.txt\" u (\$1):2 w l lw 2 lc rgb \"red\" title \"Likhtman-McLeish\"
"

printf "$plot" | gnuplot


