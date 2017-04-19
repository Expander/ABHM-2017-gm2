set terminal pdf size 10cm,8cm
set logscale x
set key box top left height 0.4
set grid
set xlabel "tan {/Symbol b}"
set ylabel "a_{/Symbol m}^{SUSY} {/Symbol \264} 10^{10}"
set format x "10^{%T}"
set output "tb_resummation.pdf"

dataTB    = "TB.dat"

plot [:] \
     dataTB u 1:($4*1e10) t 'tan {/Symbol b} resummation' \
     with lines ls 1 dt 1 lw 2, \
     dataTB u 1:(($4-$5/2)*1e10):(($4+$5/2)*1e10) t '' \
     with filledcurves ls 1 fs transparent solid 0.3, \
     dataTB u 1:($2*1e10) t 'no tan {/Symbol b} resummation' \
     with lines ls 2 dt 2 lw 2, \
