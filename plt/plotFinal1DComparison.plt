set terminal png size 1600,1200 enhanced

set output 'plots/HLLCvSLIC_test2_1D.png'

test = 'test2_1D'
solvers = "HLLC SLIC"

stats './dat/output_'.word(solvers, 1).'_'.test.'.dat' nooutput

set xlabel "x"

vars = "rho u v w p Bx By Bz e"

set size 0.9,0.9
set origin 0,0
set multiplot layout 3,3 columnsfirst scale 1.1,0.9

do for[j=3:11] {
	set ylabel word(vars, (j-2))

	plot './dat/output_'.word(solvers, 1).'_'.test.'.dat' index (11) using 1:j t \
	     word(solvers, 1) pt 7 ps 0.01, \
	     './dat/output_'.word(solvers, 2).'_'.test.'.dat' index (11) using 1:j t \
	     word(solvers, 2) pt 7 ps 0.01
}
unset multiplot

