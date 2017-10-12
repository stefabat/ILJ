set terminal cairolatex pdf size 12cm,8cm
set output 'n3_t55x_ab-initio_eint.tex'

set xrange [1:15]
set yrange [-60:-20]
set xtics 1,2,15
set xtics nomirror
set ytics nomirror
set xlabel '\# units'
set ylabel 'Interaction energy (kcal/mol)'
set key samplen 3 invert spacing 1.5
set border lw 2

# Fitted values already in from Julia
mp2_fit(x) = a*x / (b+x); a = -64.7483; b = 4.10988;
 cc_fit(x) = c*x / (d+x); c = -65.9295; d = 2.43184;

# Fitting starting with the same points as in Julia
#a = -65.0; b = 3.5; c = -65.0; d = 3.5;
#fit mp2_fit(x) 'n3_t55x_ri-scs-mp2_mix-cc-pvdz_eint.dat'  u 1:2 via a,b
#fit  cc_fit(x) 'n3_t55x_dlpno-ccsdt_eint.dat' 			  u 1:3 via c,d

#stats 'n3_t55x_ri-scs-mp2_mix-cc-pvdz_eint.dat'  u 1:2 name 'mp2'
#stats 'n3_t55x_dlpno-ccsdt_meint.dat' 			 u 1:3 name 'cc'

set label "$\\to -65.93$" at 15.1,-57

plot mp2_fit(x) lc rgb 'gray' lw 3 t '', \
	  cc_fit(x) lc rgb 'gray' lw 3 t '', \
	 'n3_t55x_ri-scs-mp2_mix-cc-pvdz_eint.dat' u 1:2 w p pt 4 lc rgb 'light-blue' lw 5 t 'RI-SCS-MP2/cc-pVDZ', \
     'n3_t55x_ri-scs-mp2_mix-cc-pvtz_eint.dat' u 1:2 w p pt 5 lc rgb 'skyblue'    lw 5 t 'RI-SCS-MP2/cc-pVTZ', \
     'n3_t55x_dlpno-ccsdt_eint.dat' 		   u 1:2 w p pt 8 lc rgb 'red'        lw 5 t 'DLPNO-CCSD(T)/cc-pVDZ', \
     'n3_t55x_dlpno-ccsdt_eint.dat' 		   u 1:3 w p pt 9 lc rgb 'dark-red'   lw 5 t 'approximate DLPNO-CCSD(T)/cc-pVTZ'

set output


