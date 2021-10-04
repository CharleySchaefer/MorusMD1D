#!/bin/bash

#=========================================
# USER INPUT
  # Program & numerical properties
dosim=1
dopp=1
verbose=1                  # Show shell output
outdir0="demo/StickyDiffusivity"
Ntime=8000000
NseedFirst=1
NseedLast=25
time_fac=0.4

  # Polymer properties
ZE=10                      # Number of entanglements per chain
# ZS=2
sequence="00000100000"               # ZS=1; N=11
#sequence="000000000010000000000"    # ZS=1; N=21
# ZS=2
#sequence="00100000100"              # ZS=2; N=11
#sequence="000001000000000100000"    # ZS=2; N=21
# ZS=5
#sequence="01010101010"              # ZS=5; N=11
#sequence="000100010010010001000"    # ZS=5; N=21
#sequence="001001001001001000100100100100100"  # used for ZS=10; NS=10; N=33; Ntime=80000000

  # sticker properties
#fraction_closed=0.9        # Fraction of closed stickers
sticker_kdiss=0.0001         # Sticker dissociation rate
sticker_ksb=0              # Rate of bondswapping

# END USER INPUT
#=========================================


for i in `seq 1 1` ; do

fraction_closed="0."$i
#=========================================
# CREATE SUB DIRECTORY
# Count number of stickers (ZS) in sequence
ZS=`awk -vSeq=$sequence 'BEGIN{ZS=0; N=length(Seq); split(Seq, chars, ""); for (i=0; i<N;i++) {if(chars[i]=='1'){ZS+=1}}; print ZS}'`

Nbeads=`awk -vSeq=$sequence 'BEGIN{ print length(Seq) }'`
subdir="ZE"$ZE"ZS"$ZS"N"$Nbeads"p"$fraction_closed"kdiss"$sticker_kdiss"ksb"$sticker_ksb"dt"$time_fac


outdir=$outdir0/$subdir
mkdir -p $demo
mkdir -p $outdir0
mkdir -p $outdir
cp $0  $outdir

if [ $verbose ] ; then 
  echo "ZS="$ZS 
  echo "Nbeads="$Nbeads 
  echo "data will be exported to sub directory "\"$subdir\" 
fi
#=========================================

if [ $dosim -eq 1  ] ; then 
  for rseed in `seq $NseedFirst $NseedLast` ; do
    mkdir -p $outdir/seed$rseed 
./MorusMD --outdir $outdir/seed$rseed \
           --rseed $rseed --time-fac $time_fac \
           --Ntime $Ntime \
           --export-logt \
           --Ze $ZE \
           --sequence $sequence \
           --sticker-p $fraction_closed \
           --sticker-kdiss $sticker_kdiss \
           --sticker-kbondswap $sticker_ksb
  done
fi

# Postprocess 
utils="../utils"
if [ $dopp -eq 1 ] ; then
echo "Average timeprogress (using .m script)"
$utils/average_timeprogress.sh $outdir
echo "Diffusivity (saved to $subdir/diffusivity.out)":

$utils/getDiffusivityFromAveragedTimeProgress.sh $outdir/timeprogress_averaged.out $ZE > $outdir/diffusivity.out

echo "  D: `awk '{print $1" +/- "$2}' "$outdir/diffusivity.out"`"

plot="set terminal pngcairo enhanced
set output \"$outdir/rmsd.png\"

set xlabel \"{/Times-Italic t}/{/Symbol t}_R\"
set log x
set ylabel \"RMSD/a\"
set log y 

PI=3.1415
DR=1.0/(3*PI**2*$ZE) # Rouse diffusivity of non-sticky chain

set nokey

plot \\
\"$outdir/timeprogress_averaged.out\" u 1:2:3 w e ps 1.2 pt 7, \\
sqrt(2*DR*x)
"
printf "$plot" > $outdir/plot_rmsd.plt
printf "$plot" | gnuplot


fi
done
