set terminal pngcairo enhanced
set output "PD.png"


# sqZi = 1/Z
Tf(sqZi) =(1+sqZi*(-2*1.69*+sqZi*(4.17-1.55*sqZi)));

Ze=10

tauRouse(Zs)=Zs**2
tauRep(Ze, Zs)=3*Tf(1/Zs)*Ze*tauRouse(Zs)

set xrange [1:30]
set log x
set xlabel "Zs"
set log y
set ylabel "~{/Symbol g}{1.1.}{/Symbol t}_S"
set key reverse left top
set samples 1000
set nokey
plot \
'+' using ($1>3.61? $1 : 1/0):(1/tauRouse($1)):(1/tauRep(10, $1)) with filledcurves closed fc rgb "#FFEEFF",\
1/tauRouse(x) title "Sticky Rouse time",\
(x>3.61?1/tauRep(10,x):1/0) title "Sticky Reptation time, Ze=10"

