set terminal pdf size 8cm,8cm
set grid
set ylabel "a_{/Symbol m}^{SUSY} {/Symbol \264} 10^{10}"

set key box top left height 0.4 width -2
set output "OS-vs-DR_splitting.pdf"
set xlabel "M_2 = m_A = m_l / GeV, m_r = 200 GeV"
data = "scan_MS_OS-vs-DR_splitting.dat"

plot [:] \
     data u 1:($4*1e10) t 'GM2Calc (OS)' with lines ls 1 dt 1 lw 2, \
     data u 1:(($4-$5/2)*1e10):(($4+$5/2)*1e10) t 'GM2Calc uncertainty' \
     with filledcurves ls 1 fs transparent solid 0.3, \
     data u 1:($2*1e10) t 'FlexibleSUSY (DR)' with lines ls 2 dt 2 lw 2, \
     data u 1:(($2-$3/2)*1e10):(($2+$3/2)*1e10) t 'scale variation' \
     with filledcurves ls 2 fs transparent solid 0.3, \

set key box bottom right height 0.4 width -2
set output "OS-vs-DR.pdf"
set xlabel "M_2 = m_A = m_l = m_r / GeV"
data = "scan_MS_OS-vs-DR.dat"

replot