#!/bin/bash
#
# Consider a tube at time 0. The fraction of the tube not relaxed at time t
# is mu(t), see Doi 1983. The software has two methods to calculate mu(t)
# during the simulation:
#
# I.  (see details Likhtman-McLeish, Macromolecules 2002): space is discretised,
#     and the statistics is collected by taking data for each 'grid segment'.
#     This gives good statistics, and is a bit more computationally demanding 
#     than method II. Method I requires input --mu-Nx and --mu-dx to discretise 
#     space. If these are not given as argument to the program, method I 
#     is not used.
# II. At time t0 the chain ends and tube length are remembered, and at time
#     t the fraction of the tube length, mu(t), that has not yet relaxed is 
#     tracked, see e.g. Doi 1983. This method is not computationally demanding,
#     but gives poorer statistics than method I.
#
# The mu(t) calculation is controlled by the following input to the program
#   --mu-time-steps  <time step interval between mu(t) calculations>
#   Time binning:
#   --mu-Nbins <n> --mu-time-L <lower> --mu-time-U <upper > 
#   Method I discretisation of space:
#   --mu-dx <dx> --mu-Nx <Nx> 

# Simulation settings
do_simulate=1
Ze=10             # Number of entanglements per chain
Nbeads=100         # Number of beads in first simulation
Niter=200000      # Number of time steps
Nprint=2000      # Print timeprogress after this number of steps
time_fac=0.4      # Time step size (should be <0.5) 

# Calculation of mu(t)
mu_N=1            # Calculate statistics after this number of steps

# mu(t) - binning of time interval
tL=1e-3        # mu(t) histogram: lower time boundary
tU=1e6        # mu(t) histogram: upper time boundary
nbins=600     # mu(t) histogram: number of time bins

# mu(t) - discretisation of space (only if method 1 is activated)
mu_dx=0.1 # units - size of entanglement blob
mu_Nx=2000 # Nx*dx should be much larger than the chain length Ze 

# Postprocessing
do_postprocess=1

#-----------------------------------------------------
# FIRST SIMULATION (with NbeadsA)
Sequence=""
for i in `seq 1 $Nbeads`; do
  Sequence=$Sequence"0"
done
outdir=demo/demo_calculate_mu

echo "Start simulation - data will be exported to $outdir "
mkdir -p demo
mkdir -p $outdir

if [ $do_simulate -eq 1 ] ; then
#valgrind  --track-origins=yes 
./MorusMD --outdir $outdir  --sequence $Sequence --Ze $Ze \
	   --Ntime $Niter \
	   --set-boundary B \
	   --export-logt --Nprint $Nprint \
	   --mu-time-steps  $mu_N \
	   --mu-Nbins $nbins --mu-time-L $tL --mu-time-U $tU \
	   --mu-dx $mu_dx --mu-Nx $mu_Nx \
	   --time-fac $time_fac > $outdir/mu.out
echo "Simulation finished."
fi


#-----------------------------------------------------
# POSTPROCESS: Compare both simulations with theory
#              conclusion: improved agreement with theory
#                          for with increased Nbeads
if [ $do_postprocess -eq 1 ] ; then
echo "Postprocessing"
utils="../utils"
echo "  Call get_mu_LMtheory.m to calculate theoretical mu(t)"
$utils/run_mscript.sh "addpath('$utils/Octave') ; get_mu_LMtheory($Ze, $tL, $tU, $nbins)" > $outdir/LMtheory.out
 
echo "  theoretical mu(t) exported to $outdir/LMtheory.out"
echo "  Call fit_dmu_LMtheory.m to extract simulated mu(t) parameters Cmu, Gf, Tf."
run_in_matlab=`$utils/is_matlab_installed.sh` # 0: matlab not installed; 1: matlab is installed
$utils/run_mscript.sh "addpath('$utils/Octave') ; fit_dmu_LMtheory('demo/demo_calculate_mu/mu.out', $Ze, $run_in_matlab)" 



echo "Plotting (gnuplot)"

plot="set terminal pngcairo enhanced
set output \"$outdir/mu.png\"

set xlabel \"t/{/Symbol t}_e\"
set log x
set xtics format \"10^{%%T}\"

set ylabel \"{/Symbol m}(t)\"
set yrange [0:1]

set key left bottom 
set grid lw 0.5 lt 1 lc rgb \"grey\"

plot \\
\"$outdir/mu.out\" u (\$1):4 w p ps 1.1 pt 6 lc rgb \"#6666DD\" title \"N=$Nbeads; Method 1\",\\
\"$outdir/mu.out\" u (\$1):2 w p ps 1.1 pt 6 lc rgb \"#DD6666\" title \"N=$Nbeads; Method 2\",\
\"$outdir/LMtheory.out\" u (\$1):2 w l lw 2 lc rgb \"black\" title \"Likhtman-McLeish\"
"

printf "$plot" | gnuplot
plot="set terminal pngcairo enhanced
set output \"$outdir/dmu.png\"

set xlabel \"t/{/Symbol t}_e\"
set log x
set xtics format \"10^{%%T}\"

set ylabel \"d{/Symbol m}(t)\"
set yrange [0:]

set key left bottom 
set grid lw 0.5 lt 1 lc rgb \"grey\"

plot \\
\"$outdir/mu.out\" u (\$1):5 w p ps 1.1 pt 6 lc rgb \"#6666DD\" title \"N=$Nbeads; Method 1\",\\
\"$outdir/mu.out\" u (\$1):3 w p ps 1.1 pt 6 lc rgb \"#DD6666\" title \"N=$Nbeads; Method 2\",\
\"$outdir/LMtheory.out\" u (\$1):3 w l lw 2 lc rgb \"black\" title \"Likhtman-McLeish\"
"

printf "$plot" | gnuplot
fi
echo "End of demo"

