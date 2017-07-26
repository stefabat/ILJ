set terminal cairolatex pdf size 9cm,6cm
set output 'n3_t55x_eint.tex'

set xrange [2:14]
set yrange [-55:-20]
set xtics 3,2,13
set xtics nomirror
set ytics nomirror
set xlabel '\# units'
set ylabel 'Interaction energy (kcal/mol)'
set key samplen 3
set border lw 2

plot 'n3_t55x_dlpno-ccsdt_mix-cc-pvdtz_eint.dat' u 1:2 w lp dt 2 pt 9 lc rgb 'skyblue' lw 5 t 'DLPNO-CCSD(T)/DZ', \
     'n3_t55x_ri-scs-mp2_mix-cc-pvdtz_eint.dat' u 1:2 w lp dt 2 pt 7 lc rgb 'dark-turquoise' lw 5 t 'SCS-MP2/DZ', \
     'n3_t55x_ri-b97d3_mix-cc-pvdtz_eint.dat' u 1:2 w lp dt 2 pt 5 lc rgb 'midnight-blue' lw 5 t 'B97D3/DZ', \
     'n3_t55x_ri-scs-mp2_mix-cc-pvtz_eint.dat' u 1:2 w lp dt 2 pt 6 lc rgb 'orange-red' lw 5 t 'SCS-MP2/TZ', \
     'n3_t55x_ri-b97d3_mix-cc-pvtz_eint.dat' u 1:2 w lp dt 2 pt 4 lc rgb 'coral' lw 5 t 'B97D3/TZ'

set output


