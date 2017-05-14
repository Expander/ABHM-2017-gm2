set terminal pdf size 10cm,8cm
set logscale x
set key box top left height 0.4
set grid
set xlabel "tan {/Symbol b}"
set ylabel "a_{/Symbol m}^{SUSY} {/Symbol \264} 10^{10}"
set format x "10^{%T}"

set style line 1 lc rgb '#ff0000' lt 1 lw 2 dt 1
set style line 2 lc rgb '#0000ff' lt 2 lw 2 dt 2
set style line 3 lc rgb '#00aa00' lt 3 lw 2 dt 4
set style line 4 lc rgb '#999999' lt 3 lw 2 dt 4

lt(x,val) = x < val ? x : 1/0

set output "tb_resummation.pdf"
dataTB    = "TB.dat"

plot [:] \
     dataTB u 1:($4*1e10) t 'tan {/Symbol b} resummation' \
     with lines ls 1 dt 1 lw 2, \
     dataTB u (lt($1,100)):($2*1e10) t 'no tan {/Symbol b} resummation' \
     with lines ls 2 dt 2 lw 2, \
     # dataTB u 1:($6*1e10) t 'FeynHiggs' \
     # with lines ls 3 dt 4 lw 2, \

set output "tb_resummation_uncertainty.pdf"

plot [:] \
     dataTB u 1:($4*1e10) t 'tan {/Symbol b} resummation' \
     with lines ls 1 dt 1 lw 2, \
     dataTB u 1:(($4-$5/2)*1e10):(($4+$5/2)*1e10) t '' \
     with filledcurves ls 1 fs transparent solid 0.3, \
     dataTB u (lt($1,100)):($2*1e10) t 'no tan {/Symbol b} resummation' \
     with lines ls 2 dt 2 lw 2, \
     # dataTB u 1:($6*1e10) t 'FeynHiggs' \
     # with lines ls 3 dt 4 lw 2, \
