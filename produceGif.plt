set terminal gif animate delay 3
set output filename.'.gif'
stats './dat/'.filename.'.dat' nooutput

do for [i=1:int(STATS_blocks)] {
	set size 1,1
	set origin 0,0
	set multiplot layout 2,2 columnsfirst scale 1.1,0.9
	    plot './dat/'.filename.'.dat' index (i-1) using 1:2 t "rho" with lines
	    plot './dat/'.filename.'.dat' index (i-1) using 1:3 t "u" with lines
	    plot './dat/'.filename.'.dat' index (i-1) using 1:4 t "p" with lines
	    plot './dat/'.filename.'.dat' index (i-1) using 1:5 t "e" with lines
	unset multiplot
}
