set terminal cairolatex pdf size 12cm,10cm
set output 'n3_benzene_pes.tex'

set xrange [3:8]
set yrange [-200:200]
#set xtics 1,2,15
#set ytics -60,10,40
set xtics nomirror
set ytics nomirror
set xlabel 'R [\AA]'
set ylabel 'Interaction energy [meV]'
set key samplen 3
set border lw 2
set grid
set xzeroaxis linetype -1


plot 'n3_benzene.dat' u 1:2 w l lc rgb 'blue' lw 3 t '$V_{ilj}$', \
     'n3_benzene.dat' u 1:3 w l lc rgb 'orange' lw 3 t '$V_{ind}$', \
     'n3_benzene.dat' u 1:4 w l lc rgb 'green' lw 3 t '$V_{els}$', \
     'n3_benzene.dat' u 1:5 w l lc rgb 'red' lw 6 t '$V_{tot}$'


