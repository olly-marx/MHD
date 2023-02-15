set terminal png size 1600,1200 enhanced
set output 'plots/'.filename.'_all.png'

stats './dat/'.filename.'.dat' nooutput

set xlabel "x"
set ylabel "y"

set pm3d
set palette rgbformulae 33,13,10

set ticslevel 0
set hidden3d

vars = "rho u v p e"

set size 1,1
set origin 0,0
set multiplot layout 3,2 columnsfirst scale 1.1,0.9
do for[j=3:7] {
		splot './dat/'.filename.'.dat' index (7) using 1:2:j \
		t word(vars, (j-2)) with pm3d
}
unset multiplot

