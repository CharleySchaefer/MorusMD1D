#!/bin/bash

# PURPOSE: 
# Simulate sticky reptation in absence of flow
# for a single sequence, and fixed physical properties
# such as fraction of closed stickers, and sticker lifetime

#=========================================
# USER INPUT
  # Program & numerical properties
dosim=1
do_postprocess=1
verbose=0                  # Show shell output
outdir0="demo/StickyDumbbell"
#Ntime=80000000
#Ntime=80000000
#Nprint=200
Ntime=200000000
Nprint=500000
rseed=1
#NseedFirst=1
#NseedLast=1
time_fac=0.0002
Wi=0.004

  # Polymer properties
ZE=10                      # Number of entanglements per chain
#ZS=2
#sequence="11" # STICKY DUMBBELL
#sequence="010010"
#sequence="01000010"
sequence="010000000010"
#sequence="01000000000000000010"
#sequence="010000000000000000000000000000000010"
#sequence="01000000000000000000000000000000000000000000000000000000000000000010"
#ZS=10
#sequence="000010000100001000010000100001000010000100001000010000"    # Polymer sequence (1=sticky bead, 0=non-sticky bead)

  # sticker properties
fraction_closed=0.9       # Fraction of closed stickers
sticker_kdiss=0.001         # Sticker dissociation rate
sticker_ksb=0              # Rate of bondswapping

# END USER INPUT
#=========================================
# 0.001 0.002 0.003 0.004 0.005 0.01 0.02 0.05 
for Wi in 0.03  ; do
  echo "Wi=$Wi"


#=========================================
# CREATE SUB DIRECTORY
# Count number of stickers (ZS) in sequence
ZS=`awk -vSeq=$sequence 'BEGIN{ZS=0; N=length(Seq); split(Seq, chars, ""); for (i=0; i<=N;i++) {if(chars[i]=='1'){ZS+=1}}; print ZS}'`

Nbeads=`awk -vSeq=$sequence 'BEGIN{ print length(Seq) }'`
subdir="ZE"$ZE"ZS"$ZS"N"$Nbeads"p"$fraction_closed"kdiss"$sticker_kdiss"ksb"$sticker_ksb"dt"$time_fac"Wi"$Wi



outdir=$outdir0/$subdir
mkdir -p "demo"
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
  #for rseed in `seq $NseedFirst $NseedLast` ; do
    mkdir -p $outdir 
#           --set-boundary A \
./MorusMD --outdir $outdir \
           --rseed $rseed --time-fac $time_fac \
           --Ntime $Ntime \
           --Nprint $Nprint \
           --Ze $ZE \
           --sequence $sequence \
           --sticker-p $fraction_closed \
           --sticker-kdiss $sticker_kdiss \
           --sticker-kbondswap $sticker_ksb \
           --Wi $Wi \
           --sticker-sync #\
#           --exclude-fluctuations
#           --export-logt \
#           --export-sequence \
#  done
fi

# Postprocess
if [ $do_postprocess -eq 1  ] ; then

  # GET STRETCH AND STRESS DISTRIBUTIONS
  ../utils/analyseStretch.sh $outdir
  ../utils/analyseStress.sh $outdir

  # Get sticky diffusion coefficient (units of a^2/taue)
  tauS=`awk -vkopen=$sticker_kdiss 'BEGIN {print 1.0/kopen}' `
  DSR=`../utils/get_StickyDiffusionCoefficient.m $ZE $ZS $fraction_closed $tauS`
  tauSR=`awk -vDSR=$DSR 'BEGIN {print 1.0/DSR}' `

  ExtensionRate=`awk -vtauSR=$tauSR -vWi=$Wi -vZE=$ZE 'BEGIN {print Wi*tauSR/ZE**2}'` # units of tauSR


echo "Sticker lifetime:  $tauS / units of tauE"
echo "Sticky Rouse time: $tauSR / units of tauE"
echo "Extension Rate:    $ExtensionRate"
# Plot results
plot="set terminal pngcairo enhanced dashed
set output \"$outdir/ChainCoordinates.png\"

set log x
set xlabel \"time/{/Symbol t}_{SR}\"
set xtics format \"10^{%%T}\" 
#set xrange [0: 120]

set ylabel \"x/a\"

Wi=$Wi
Z=$ZE
pi=3.1415

if (Wi==0) { \\
  Zstar=Z/2\\
} else {\\
  Zstar=sin( pi*sqrt((Wi/Z**2))*(  Z/2  )  )/(pi*sqrt((Wi/Z**2))*cos( pi*sqrt((Wi/Z**2))*(Z/2))) }
set grid lt 1 lw 0.5 lc rgb \"grey\" 

set style line 1 lt 1 lw 3 lc rgb \"black\"
set style line 2 lt 1 lw 1 lc rgb \"grey\"
set style line 3 lt 1 lw 3 lc rgb \"red\"
set key left bottom reverse

plot \"$outdir/timeprogress.out\" u (\$0/0):(\$0/0)"

# PLOT BEADS
for n in `seq 1 $Nbeads` ; do
  ycol=$(( n + 5 ))
plot="$plot,\\
\"\" u (\$1/$tauSR):(\$$ycol )  w l ls 2 notitle "
done


# HIGHLIGHT STICKERS
for n in `seq 1 $Nbeads` ; do
  ycol=$(( n + 5 ))

  issticker=${sequence:$n:1}
  if [[ $issticker == "1" ]] ; then 
plot="$plot,\\
\"\" u  (\$1/$tauSR):(\$$ycol)  w l ls 1 notitle "
  fi
done

# HIGHLIGHT END GROUPS
#  END GROUPS
plot="$plot,\\
\"$outdir/timeprogress.out\" u  (\$1/$tauSR):(\$2) w l ls 3 title \"end bead\", \\
\"\" u  (\$1/$tauSR):(\$3)  w l ls 3 notitle"

# PLOT EXPECTATION VALUE
plot="$plot, \\
 Zstar w l lw 2 lt 2 lc rgb \"blue\" title \"analytical steady state\", \\
-Zstar w l lw 2 lt 2 lc rgb \"blue\" notitle,\\
Zstar+sqrt(Zstar/3) w l lw 1 lt 2 lc rgb \"blue\" notitle,\\
Zstar-sqrt(Zstar/3) w l lw 1 lt 2 lc rgb \"blue\" notitle,\\
-Zstar+sqrt(Zstar/3) w l lw 1 lt 2 lc rgb \"blue\" notitle,\\
-Zstar-sqrt(Zstar/3) w l lw 1 lt 2 lc rgb \"blue\" notitle
"
printf "$plot" | gnuplot

# Plot stretch transient
fac=`awk -veps=$ExtensionRate -vp=$fraction_closed 'BEGIN {printf eps*p}'`
plot="set terminal pngcairo enhanced
set output \"$outdir/StretchTransient.png\"

set style line 1 lt 1 lw 3 lc rgb \"black\"
set style line 2 lt 1 lw 1 lc rgb \"grey\"
set style line 3 lt 1 lw 3 pt 6 lc rgb \"red\"

set xlabel \"time/{/Symbol t}_{SR}\"

set ylabel \"{/Symbol l}\"
set log y
set ytics format \"10^{%%T}\"

set grid lt 1 lw 0.5 lc rgb \"grey\" 

set key left top

plot \"$outdir/timeprogress.out\" u (\$1/$tauSR):( (\$3-\$2)/$ZE) w p ls 3 title \"~{/Symbol e}{1.1.}p{/Symbol t}_{SR}=$fac\"#,\\
#exp( $fraction_closed*$ExtensionRate*x/($ZS+1) ) w l ls 1 title \"exp( ~{/Symbol e}{1.1.} {/Symbol t}_{SR} p/(Z_{s}+1) )\"

"
printf "$plot" | gnuplot

# Plot stress transient
fac=`awk -veps=$ExtensionRate -vp=$fraction_closed 'BEGIN {printf eps*p}'`
plot="set terminal pngcairo enhanced
set output \"$outdir/StressTransient.png\"

set style line 1 lt 1 lw 3 lc rgb \"black\"
set style line 2 lt 1 lw 1 lc rgb \"grey\"
set style line 3 lt 1 lw 3 pt 6 lc rgb \"red\"

set xlabel \"time/{/Symbol t}_{SR}\"

set ylabel \"{/Symbol s}\"
set log y
set ytics format \"10^{%%T}\"

set grid lt 1 lw 0.5 lc rgb \"grey\" 

set key left top

plot \"$outdir/timeprogress.out\" u (\$1/$tauSR):( (\$5) ) w p ls 3 title \"~{/Symbol e}{1.1.}p{/Symbol t}_{SR}=$fac\"#,\\
#exp( $fraction_closed*$ExtensionRate*x/($ZS+1) ) w l ls 1 title \"exp( ~{/Symbol e}{1.1.} {/Symbol t}_{SR} p/(Z_{s}+1) )\"

"
printf "$plot" | gnuplot


# Plot stretch distribution
plot="set terminal pngcairo enhanced
set output \"$outdir/StretchDistribution.png\"

set style line 1 lt 1 lw 3 lc rgb \"black\"
set style line 2 lt 1 lw 1 lc rgb \"grey\"
set style line 3 lt 1 lw 3 pt 6 lc rgb \"red\"

set xlabel \"{/Symbol l}\"
set log x

set ylabel \"P({/Symbol l})\"
set log y
set ytics format \"10^{%%T}\"

set grid lt 1 lw 0.5 lc rgb \"grey\" 

set key left top

plot \"$outdir/StretchDistribution.out\" u 1:2 w p ls 3 notitle \"~{/Symbol e}{1.1.}p{/Symbol t}_{SR}=$fac\"#,\\
#exp( $fraction_closed*$ExtensionRate*x/($ZS+1) ) w l ls 1 title \"exp( ~{/Symbol e}{1.1.} {/Symbol t}_{SR} p/(Z_{s}+1) )\"

"
printf "$plot" | gnuplot


# Plot stress distribution
plot="set terminal pngcairo enhanced
set output \"$outdir/StressDistribution.png\"

set style line 1 lt 1 lw 3 lc rgb \"black\"
set style line 2 lt 1 lw 1 lc rgb \"grey\"
set style line 3 lt 1 lw 3 pt 6 lc rgb \"red\"

set xlabel \"{/Symbol s}\"
set log x

set ylabel \"P({/Symbol s})\"
set log y
set ytics format \"10^{%%T}\"

set grid lt 1 lw 0.5 lc rgb \"grey\" 

set key left top

plot \"$outdir/StressDistribution.out\" u 1:2 w p ls 3 notitle \"~{/Symbol e}{1.1.}p{/Symbol t}_{SR}=$fac\"#,\\
#exp( $fraction_closed*$ExtensionRate*x/($ZS+1) ) w l ls 1 title \"exp( ~{/Symbol e}{1.1.} {/Symbol t}_{SR} p/(Z_{s}+1) )\"

"
printf "$plot" | gnuplot

fi

done
