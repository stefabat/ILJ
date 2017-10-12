set terminal cairolatex pdf size 12cm,8cm
set output 'n3_t55x_comp_eint.tex'

set xrange [3:13]
set yrange [-65:-15]
set xtics 1,2,15
set xtics nomirror
set ytics nomirror
set xlabel '\# units'
set ylabel 'Interaction energy (kcal/mol)'
set key samplen 3 spacing 1.5
set border lw 2
set grid

cc_fit(x) = c*x / (d+x); c = -65.9295; d = 2.43184;

plot cc_fit(x) lc rgb 'gray' lw 6 t '\emph{ab-initio}', \
	 'n3_t55x_classical_eint.dat'  u 1:2 w l lc rgb 'red'    lw 6 t 'classical'

set output


