set terminal png size 1600,1200 enhanced
set output 'plots/'.filename.'_all.png'

stats './dat/'.filename.'.dat' nooutput

set xlabel "x"
set ylabel "y"

vars = "rho u v w pT Bx By Bz e E"

set size 1,1
set origin 0,0
set multiplot layout 3,3 columnsfirst scale 1.1,0.9
do for[j=3:11] {
		plot './dat/'.filename.'.dat' index (11) using 1:j \
		t word(vars, (j-2)) pt 7 ps 0.01
}
unset multiplot

