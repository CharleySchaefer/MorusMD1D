set terminal pngcairo enhanced dashed
set output "demo_plot.png"

set xlabel "t/{/Symbol t}_e"
set log x
set xtics format "10^{%T}"

set ylabel "d{/Symbol m}, dR"
set log y
set ytics format "10^{%T}"
set yrange [0.1:4]

set grid lw 0.5 lt 1 lc rgb "black"
set key left bottom


plot "demo_data.out" u 1:3 w l lt 2 lw 4 lc rgb "#DD4444" title "d{/Symbol m}(t)", \
"" u 1:5 w l lt 1 lw 4 lc rgb "#4444DD" title "dR(t)"


