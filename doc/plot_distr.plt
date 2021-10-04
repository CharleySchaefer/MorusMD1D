set terminal pngcairo enhanced
set output "distr.png"

a=1;


var_l(l_0,a) =  a*abs(l_0)*1.0/3

P(x, l_0, a) = \
1.0/sqrt( 2*3.14*var_l(l_0,a) )*exp( -(x-l_0)**2/(2*var_l(l_0,a)) )

#set xrange [-50:50]

set xlabel "distance from sticker, x/l_0({/Symbol e})"
set ylabel "sticker distance correlation function, g(x)"
set samples 1000
l_0=0.05
plot \
'+' using ($1):(sum [n=1:20] (P($1,  n, 1.0) + P($1,  -n, 1.0))) w l lw 3 title "stretch ratio, {/Symbol l}=1",\
'+' using ($1):(sum [n=1:20] (P($1,  n, 0.5) + P($1,  -n, 0.5))) w l title "{/Symbol l}=2",\
'+' using ($1):(sum [n=1:20] (P($1,  n, 0.2) + P($1,  -n, 0.2))) w l title "{/Symbol l}=5",\
'+' using ($1):(sum [n=1:20] (P($1,  n, 0.1) + P($1,  -n, 0.1))) w l title "{/Symbol l}=10",\
'+' using ($1):(sum [n=1:20] (P($1,  n, 0.05) + P($1,  -n, 0.05))) w l title "{/Symbol l}=20",\
'+' using ($1):(sum [n=1:20] (P($1,  n, 0.02) + P($1,  -n, 0.02))) w l title "{/Symbol l}=50"

