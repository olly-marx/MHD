set terminal gif animate size 800,600 delay 50
set output 'plots/'.filename.'_all.gif'

stats './dat/'.filename.'.dat' nooutput

set xlabel "x"
set ylabel "y"

set pm3d
set palette rgbformulae 33,13,10

set ticslevel 0
set dgrid3d 200,200
set hidden3d

do for [i=1:int(STATS_blocks)] {
	set size 1,1
	set origin 0,0
	set multiplot layout 3,2 columnsfirst scale 1.1,0.9
	    splot './dat/'.filename.'.dat' index (i-1) using 1:2:3 t "rho" with lines
	    splot './dat/'.filename.'.dat' index (i-1) using 1:2:4 t "u" with lines
	    splot './dat/'.filename.'.dat' index (i-1) using 1:2:5 t "v" with lines
	    splot './dat/'.filename.'.dat' index (i-1) using 1:2:6 t "p" with lines
	    splot './dat/'.filename.'.dat' index (i-1) using 1:2:7 t "e" with lines
	unset multiplot
}
