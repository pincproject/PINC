set termoption dash
#set title "Potential of MMO versus time" font "latin modern, 12"
set autoscale
set xlabel "Timestep (1e-8 s)" font "latin modern, 12"
set ylabel "Electric potential (V)" font "latin modern, 12"
set grid
set nokey
set style line 1 lt 1 lc rgb "black" lw 1.5
set style line 2 lt 3 lc rgb "black" lw 1.5
set style line 3 lt 5 lc rgb "black" lw 1.5
plot "pot.dat" with lines linestyle 1 title "PINC potential"


set multiplot

set origin 0,0
set size 1,1
set xtics font ", 10"
set ytics font ", 10"
plot "pot.dat" with lines linestyle 2



unset title
unset xlabel
unset ylabel
unset grid
set origin 0.4, 0.5
set size 0.45, 0.35
set xrange [30000:40000]
#set yrange [104:109]
#set yrange [99:104]
#set yrange [77:82]
#set yrange [74:79]
set yrange [-230:-200]
set xtics font ", 8"
set ytics 10 font ", 8"
replot

unset multiplot