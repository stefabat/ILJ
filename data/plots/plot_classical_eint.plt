set terminal cairolatex pdf size 12cm,12cm
set output 'n3_t55x_classical_eint.tex'

set xrange [0:13]
set yrange [-65:45]
set xtics 1,2,15
set ytics -60,10,40
set xtics nomirror
set ytics nomirror
set xlabel '\# units'
set ylabel 'Interaction energy (kcal/mol)'
set key samplen 3 above center horizontal
set border lw 2
set grid
set xzeroaxis linetype -1

cc_fit(x) = c*x / (d+x); c = -65.9295; d = 2.43184;

plot cc_fit(x) lc rgb 'gray' lw 6 t '\emph{ab-initio}', \
	 'n3_t55x_classical_eint.dat'  u 1:5 w l lc rgb 'blue'   lw 2 t 'Induction', \
	 'n3_t55x_classical_eint.dat'  u 1:4 w l lc rgb 'green'  lw 2 t 'Electrostatic', \
	 'n3_t55x_classical_eint.dat'  u 1:3 w l lc rgb 'yellow' lw 2 t 'Dispersion', \
	 'n3_t55x_classical_eint.dat'  u 1:2 w l lc rgb 'red'    lw 6 t 'Total classical'

set output


