set terminal cairolatex pdf size 12cm,8cm
set output 'n3_t55x_eint.tex'

set xrange [2:14]
set yrange [-55:-20]
set xtics 3,2,13
set xtics nomirror
set ytics nomirror
set xlabel '\# units'
set ylabel 'Interaction energy (kcal/mol)'
set key samplen 3 invert
set border lw 2

plot 'n3_t55x_ri-b97d3_eint.dat'                 u 1:2 w lp dt 2 pt 7 lc rgb 'turquoise'      lw 5 t 'B97D3/DZ', \
     'n3_t55x_ri-b97d3_eint.dat'                 u 1:3 w lp dt 2 pt 5 lc rgb 'turquoise'      lw 5 t 'B97D3/TZ', \
     'n3_t55x_ri-scs-mp2_mix-cc-pvdtz_eint.dat'  u 1:2 w lp dt 2 pt 7 lc rgb 'coral'          lw 5 t 'SCS-MP2/DZ', \
     'n3_t55x_ri-scs-mp2_mix-cc-pvtz_eint.dat'   u 1:2 w lp dt 2 pt 5 lc rgb 'coral'          lw 5 t 'SCS-MP2/TZ', \
     'n3_t55x_cbs_eint.dat'                      u 1:2 w lp dt 2 pt 9 lc rgb 'coral'          lw 5 t 'SCS-MP2/CBS', \
     'n3_t55x_dlpno-ccsdt_mix-cc-pvdtz_eint.dat' u 1:2 w lp dt 2 pt 7 lc rgb 'forest-green'   lw 5 t 'DLPNO-CCSD(T)/DZ', \
     'n3_t55x_cbs_eint.dat'                      u 1:3 w lp dt 2 pt 9 lc rgb 'forest-green'   lw 5 t 'DLPNO-CCSD(T)/CBS'

set output


