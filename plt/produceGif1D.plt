set terminal gif animate size 1600,1200 delay 50
set output 'plots/'.filename.'_all.gif'

stats './dat/'.filename.'.dat' nooutput

set xlabel "x"
set ylabel "y"

vars = "rho u v w p Bx By Bz e"

do for [i=1:int(STATS_blocks)] {
	set size 1,1
	set origin 0,0
	set multiplot layout 3,3 columnsfirst scale 1.1,0.9
	do for[j=3:11] {
			plot './dat/'.filename.'.dat' index (i-1) using 1:j \
			t word(vars, (j-2)) with lines
	}
	unset multiplot
}
