set terminal png size 1600,1200 enhanced

set output 'plots/HLLCvSLIC_test2_1D.png'

test = "test2_1D"
solvers = "HLLC SLIC"

stats './dat/output_HLLC_test2_1D.dat' nooutput

set xlabel "x"
set ylabel "y"

vars = "rho u v w p Bx By Bz e"

set size 1,1
set origin 0,0
set multiplot layout 3,3 columnsfirst scale 1.1,0.9
do for[j=3:11] {
		plot './dat/output_HLLC_test2_1D.dat' index (11) using 1:j t \
		     word(vars, (j-2)) pt 7 ps 0.01, \
		     './dat/output_SLIC_test2_1D.dat' index (11) using 1:j t \
		     word(vars, (j-2)) pt 7 ps 0.01
}
unset multiplot

