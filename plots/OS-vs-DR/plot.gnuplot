set terminal pdf size 8cm,8cm
set grid
set ylabel "a_{/Symbol m}^{SUSY} {/Symbol \264} 10^{10}"

set key box top left height 0.4 width -2
set output "OS-vs-DR_splitting.pdf"
set xlabel "M_2 = m_A = m_l / GeV, m_r = 200 GeV"
data = "scan_MS_OS-vs-DR_splitting.dat"

set style line 1 lc rgb '#ff0000' lt 1 lw 2 dt 1
set style line 2 lc rgb '#00aa00' lt 2 lw 2 dt 2
set style line 3 lc rgb '#000000' lt 3 lw 3 dt 3
set style line 4 lc rgb '#999999' lt 3 lw 3 dt 3
set style line 5 lc rgb '#0000ff' lt 2 lw 2 dt 2

plot [:] \
     data u 1:($4*1e10) t 'GM2Calc' with lines ls 1, \
     data u 1:(($4-$5/2)*1e10):(($4+$5/2)*1e10) t 'GM2Calc uncertainty' \
     with filledcurves ls 1 fs transparent solid 0.3, \
     data u 1:($2*1e10) t 'FlexibleSUSY' with lines ls 2, \
     data u 1:(($2-$3/2)*1e10):(($2+$3/2)*1e10) t 'scale variation' \
     with filledcurves ls 2 lw 0 dt 1 fs transparent solid 0.3, \
     data u 1:($7*1e10) t 'FeynHiggs (OS)' with lines ls 5, \
     data u 1:($6*1e10) t 'FeynHiggs (Q = M_S)' with lines ls 3, \
     data u 1:($8*1e10) t 'FeynHiggs (Q = 455 GeV)' with lines ls 4, \

set key box bottom right height 0.4 width -2
set output "OS-vs-DR.pdf"
set xlabel "M_2 = m_A = m_l = m_r / GeV"
data = "scan_MS_OS-vs-DR.dat"

plot [:] \
     data u 1:($4*1e10) t 'GM2Calc' with lines ls 1, \
     data u 1:(($4-$5/2)*1e10):(($4+$5/2)*1e10) t 'GM2Calc uncertainty' \
     with filledcurves ls 1 fs transparent solid 0.3, \
     data u 1:($2*1e10) t 'FlexibleSUSY' with lines ls 2, \
     data u 1:(($2-$3/2)*1e10):(($2+$3/2)*1e10) t 'scale variation' \
     with filledcurves ls 2 lw 0 dt 1 fs transparent solid 0.3, \
     data u 1:($7*1e10) t 'FeynHiggs (OS)' with lines ls 5, \

# replot
