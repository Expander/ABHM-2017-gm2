set terminal pdf size 8cm,8cm
set grid
set ylabel "a_{/Symbol m}^{SUSY} {/Symbol \264} 10^{10}"

set key box top left height 0.4 width -2
set output "OS-vs-DR_splitting.pdf"
set xlabel "M_2 = m_A = m_l / GeV, m_r = 200 GeV"
data = "scan_MS_OS-vs-DR_splitting.dat"

set style line 1 lc rgb '#ff0000' lt 1 lw 2 dt 1
set style line 2 lc rgb '#00aa00' lt 2 lw 2 dt 2
set style line 3 lc rgb '#0000ff' lt 3 lw 2 dt 4

plot [:] \
     data u 1:($4*1e10) t 'GM2Calc' with lines ls 1, \
     data u 1:(($4-$5/2)*1e10):(($4+$5/2)*1e10) t 'GM2Calc uncertainty' \
     with filledcurves ls 1 fs transparent solid 0.3, \
     data u 1:($2*1e10) t 'FlexibleSUSY' with lines ls 2, \
     data u 1:(($2-$3/2)*1e10):(($2+$3/2)*1e10) t 'scale variation' \
     with filledcurves ls 2 fs transparent solid 0.3, \
     data u 1:($6*1e10) t 'FeynHiggs' with lines ls 3, \

set key box bottom right height 0.4 width -2
set output "OS-vs-DR.pdf"
set xlabel "M_2 = m_A = m_l = m_r / GeV"
data = "scan_MS_OS-vs-DR.dat"

replot
