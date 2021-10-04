set terminal pngcairo enhanced
set output "TimeScales.png"


# sqZi = 1/Z
Tf(sqZi) =(1+sqZi*(-2*1.69*+sqZi*(4.17-1.55*sqZi)));

Ze=10

tauRouse(Zs)=Zs**2
tauRep(Ze, Zs)=3*Tf(1/Zs)*Ze*tauRouse(Zs)

set xrange [1:100]
set log x
set xlabel "Zs"
set log y
set ylabel "Time Scale"
set key reverse left top
plot \
tauRouse(x) title "Sticky Rouse time",\
tauRep(2,x) title "Sticky Reptation time, Ze=2",\
tauRep(5,x) title "Sticky Reptation time, Ze=5",\
tauRep(10,x) title "Sticky Reptation time, Ze=10",\
tauRep(20,x) title "Sticky Reptation time, Ze=20",\
tauRep(50,x) title "Sticky Reptation time, Ze=50"

