set terminal postscript portrait enhanced color
set fit quiet
set size 3.0/3.0, 4.0/4.0
set size ratio 0.3/1
set tics scale 0.5

set output 'results/1D_Brio&Wu.eps'

set xtics font ",10"
set ytics font ",10"

set lmargin 10
set bmargin 3

#--------------------------------------------------------------------------------
# Print the first test Brio & Wu 1D

set multiplot layout 5,1 rowsfirst

set ylabel '{/Symbol r}' font ",16"
set yrange [0.0:1.1]
plot './dat/output_1D_Brio&Wu.dat' index (0) u 1:3 t "" w lines lc rgb \
	"#000000", \
     './dat/output_1D_Brio&Wu.dat' index (1) u 1:3 t "" w points pt 5 \
	ps 0.1

set ylabel '{/:Italic u}' font ",16"
set yrange [-0.3:0.7]
plot './dat/output_1D_Brio&Wu.dat' index (0) u 1:4 t "" w lines \
	 lc rgb "#000000", \
     './dat/output_1D_Brio&Wu.dat' index (1) u 1:4 t "" w points pt 5 \
	ps 0.1

set ylabel '{/:Italic v}' font ",16"
set yrange [-1.8:0.2]
plot './dat/output_1D_Brio&Wu.dat' index (0) u 1:5 t "" w lines \
	lc rgb "#000000", \
     './dat/output_1D_Brio&Wu.dat' index (1) u 1:5 t "" w points pt 5 \
	ps 0.1

set ylabel '{/:Italic p_T}' font ",16"
set yrange [0.6:2.0]
plot './dat/output_1D_Brio&Wu.dat' index (0) u 1:7 t "" w lines \
	lc rgb "#000000", \
     './dat/output_1D_Brio&Wu.dat' index (1) u 1:7 t "" w points pt 5 \
	ps 0.1

set xlabel '{/:Italic x}'
set ylabel '{/:Italic B_y}' font ",16"
set yrange [-1.25:1.25]
plot './dat/output_1D_Brio&Wu.dat' index (0) u 1:9 t "" w lines \
	lc rgb "#000000", \
     './dat/output_1D_Brio&Wu.dat' index (1) u 1:9 t "" w points pt 5 \
	ps 0.1

unset multiplot
unset xlabel

#--------------------------------------------------------------------------------
# Print the second test Sod Euler x=.5

set size ratio 0.55/1

set output 'results/2D_Sod_Euler_x.eps'

set multiplot layout 4,1 rowsfirst

set ylabel '{/:Italic {/Symbol r}}' font ",16"
set yrange [0.0:1.1]
plot './dat/output_2D_Sod_Euler_x=.5.dat' index (0) u 1:3 t "" w lines lc rgb \
	"#000000", \
     './dat/output_2D_Sod_Euler_x=.5.dat' index (1) u 1:3 t "" w points pt 5 \
	ps 0.1

set ylabel '{/:Italic u}' font ",16"
set yrange [-0.1:1.0]
plot './dat/output_2D_Sod_Euler_x=.5.dat' index (0) u 1:4 t "" w lines \
	 lc rgb "#000000", \
     './dat/output_2D_Sod_Euler_x=.5.dat' index (1) u 1:4 t "" w points pt 5 \
	ps 0.1

set ylabel 'p' font ",16"
set yrange [0.0:1.1]
plot './dat/output_2D_Sod_Euler_x=.5.dat' index (0) u 1:7 t "" w lines \
	lc rgb "#000000", \
     './dat/output_2D_Sod_Euler_x=.5.dat' index (1) u 1:7 t "" w points pt 5 \
	ps 0.1

set xlabel '{/:Italic x}' font ",16"
set ylabel '{/:Italic {/Symbol e}}' font ",16"
set yrange [1.7:2.95]
plot './dat/output_2D_Sod_Euler_x=.5.dat' index (0) u 1:11 t "" w lines \
	lc rgb "#000000", \
     './dat/output_2D_Sod_Euler_x=.5.dat' index (1) u 1:11 t "" w points pt 5 \
	ps 0.1

unset multiplot
unset xlabel

#--------------------------------------------------------------------------------
# Print the third test with y=.5

set size ratio 0.55/1

set output 'results/2D_Sod_Euler_y.eps'

set multiplot layout 4,1 rowsfirst

set ylabel '{/:Italic {/Symbol r}}' font ",16"
set yrange [0.0:1.1]
plot './dat/output_2D_Sod_Euler_y=.5.dat' index (0) u 2:3 t "" w lines lc rgb \
	"#000000", \
     './dat/output_2D_Sod_Euler_y=.5.dat' index (1) u 2:3 t "" w points pt 5 \
	ps 0.1

set ylabel '{/:Italic v}' font ",16"
set yrange [-0.1:1.0]
plot './dat/output_2D_Sod_Euler_y=.5.dat' index (0) u 2:5 t "" w lines \
	 lc rgb "#000000", \
     './dat/output_2D_Sod_Euler_y=.5.dat' index (1) u 2:5 t "" w points pt 5 \
	ps 0.1

set ylabel 'p' font ",16"
set yrange [0.0:1.1]
plot './dat/output_2D_Sod_Euler_y=.5.dat' index (0) u 2:7 t "" w lines \
	lc rgb "#000000", \
     './dat/output_2D_Sod_Euler_y=.5.dat' index (1) u 2:7 t "" w points pt 5 \
	ps 0.1

set xlabel '{/:Italic y}' font ",16"
set ylabel '{/:Italic {/Symbol e}}' font ",16"
set yrange [1.7:2.95]
plot './dat/output_2D_Sod_Euler_y=.5.dat' index (0) u 2:11 t "" w lines \
	lc rgb "#000000", \
     './dat/output_2D_Sod_Euler_y=.5.dat' index (1) u 2:11 t "" w points pt 5 \
	ps 0.1

unset multiplot
unset xlabel

#--------------------------------------------------------------------------------
# Print the 4th Test Sod with diagonal y=1-x discontinuity

set size ratio 0.55/1

set output 'results/2D_Sod_Euler_y=1-x.eps'

set multiplot layout 4,1 rowsfirst

set ylabel '{/:Italic {/Symbol r}}' font ",16"
set yrange [0.0:1.1]
plot './dat/output_2D_Sod_Euler_y=x.dat' index (0) u 1:3 t "" w lines lc rgb \
	"#000000", \
     './dat/output_2D_Sod_Euler_y=x.dat' index (1) u 1:3 t "" w points pt 5 \
	ps 0.1

set ylabel '||{/:Bold v}||' font ",16"
set yrange [-0.1:1.0]
plot './dat/output_2D_Sod_Euler_y=x.dat' index (0) u 1:12 t "" w lines \
	 lc rgb "#000000", \
     './dat/output_2D_Sod_Euler_y=x.dat' index (1) u 1:12 t "" w points pt 5 \
	ps 0.1

set ylabel 'p' font ",16"
set yrange [0.0:1.1]
plot './dat/output_2D_Sod_Euler_y=x.dat' index (0) u 1:7 t "" w lines \
	lc rgb "#000000", \
     './dat/output_2D_Sod_Euler_y=x.dat' index (1) u 1:7 t "" w points pt 5 \
	ps 0.1

set xlabel '{/:Italic x = y}' font ",16"
set ylabel '{/:Italic {/Symbol e}}' font ",16"
set yrange [1.7:2.95]
plot './dat/output_2D_Sod_Euler_y=x.dat' index (0) u 1:11 t "" w lines \
	lc rgb "#000000", \
     './dat/output_2D_Sod_Euler_y=x.dat' index (1) u 1:11 t "" w points pt 5 \
	ps 0.1

unset multiplot
unset xlabel

#--------------------------------------------------------------------------------
# Print the 6th Test B&W with x=.5 discontinuity

set size ratio 0.3/1

set output 'results/2D_B&W_MHD_x=.5.eps'

unset yrange
set multiplot layout 5,1 rowsfirst

set ylabel '{/Symbol r}' font ",16"
set yrange [0.0:1.1]
plot './dat/output_2D_B&W_MHD_x=.5.dat' index (0) u 1:3 t "" w lines lc rgb \
	"#000000", \
     './dat/output_2D_B&W_MHD_x=.5.dat' index (1) u 1:3 t "" w points pt 5 \
	ps 0.1

set ylabel '{/:Italic u}' font ",16"
set yrange [-0.4:0.7]
plot './dat/output_2D_B&W_MHD_x=.5.dat' index (0) u 1:4 t "" w lines \
	 lc rgb "#000000", \
     './dat/output_2D_B&W_MHD_x=.5.dat' index (1) u 1:4 t "" w points pt 5 \
	ps 0.1

set ylabel '{/:Italic v}' font ",16"
set yrange [-1.8:0.2]
plot './dat/output_2D_B&W_MHD_x=.5.dat' index (0) u 1:5 t "" w lines \
	lc rgb "#000000", \
     './dat/output_2D_B&W_MHD_x=.5.dat' index (1) u 1:5 t "" w points pt 5 \
	ps 0.1

set ylabel '{/:Italic p_T}' font ",16"
set yrange [0.6:1.9]
plot './dat/output_2D_B&W_MHD_x=.5.dat' index (0) u 1:7 t "" w lines \
	lc rgb "#000000", \
     './dat/output_2D_B&W_MHD_x=.5.dat' index (1) u 1:7 t "" w points pt 5 \
	ps 0.1

set xlabel '{/:Italic x}'
set ylabel '{/:Italic B_y}' font ",16"
set yrange [-1.15:1.15]
plot './dat/output_2D_B&W_MHD_x=.5.dat' index (0) u 1:9 t "" w lines \
	lc rgb "#000000", \
     './dat/output_2D_B&W_MHD_x=.5.dat' index (1) u 1:9 t "" w points pt 5 \
	ps 0.1

unset multiplot
unset xlabel

#--------------------------------------------------------------------------------
# Print the 7th Test B&W with y=.5 discontinuity

set size ratio 0.3/1

set output 'results/2D_B&W_MHD_y=.5.eps'

unset yrange
set multiplot layout 5,1 rowsfirst

set ylabel '{/Symbol r}' font ",16"
set yrange [0.0:1.1]
plot './dat/output_2D_B&W_MHD_y=.5.dat' index (0) u 2:3 t "" w lines lc rgb \
	"#000000", \
     './dat/output_2D_B&W_MHD_y=.5.dat' index (1) u 2:3 t "" w points pt 5 \
	ps 0.1

set ylabel '{/:Italic u}' font ",16"
set yrange [-1.8:0.2]
plot './dat/output_2D_B&W_MHD_y=.5.dat' index (0) u 2:4 t "" w lines \
	 lc rgb "#000000", \
     './dat/output_2D_B&W_MHD_y=.5.dat' index (1) u 2:4 t "" w points pt 5 \
	ps 0.1

set ylabel '{/:Italic v}' font ",16"
set yrange [-0.4:0.7]
plot './dat/output_2D_B&W_MHD_y=.5.dat' index (0) u 2:5 t "" w lines \
	lc rgb "#000000", \
     './dat/output_2D_B&W_MHD_y=.5.dat' index (1) u 2:5 t "" w points pt 5 \
	ps 0.1

set ylabel '{/:Italic p_T}' font ",16"
set yrange [0.6:1.9]
plot './dat/output_2D_B&W_MHD_y=.5.dat' index (0) u 2:7 t "" w lines \
	lc rgb "#000000", \
     './dat/output_2D_B&W_MHD_y=.5.dat' index (1) u 2:7 t "" w points pt 5 \
	ps 0.1

set xlabel '{/:Italic y}'
set ylabel '{/:Italic B_x}' font ",16"
set yrange [-1.15:1.15]
plot './dat/output_2D_B&W_MHD_y=.5.dat' index (0) u 2:8 t "" w lines \
	lc rgb "#000000", \
     './dat/output_2D_B&W_MHD_y=.5.dat' index (1) u 2:8 t "" w points pt 5 \
	ps 0.1

unset multiplot
unset xlabel

#--------------------------------------------------------------------------------
# Print the 8th test Brio and Wu with diagonal y=1-x discontinuity

set size ratio 0.3/1

set output 'results/2D_B&W_MHD_y=1-x.eps'

set multiplot layout 5,1 rowsfirst
unset yrange
set xrange [200:1000]

set ylabel '{/:Italic {/Symbol r}}' font ",16"
set yrange [0.0:1.1]
plot './dat/output_2D_B&W_MHD_y=x.dat' index (0) u 15:3 t "" w lines lc rgb \
	"#000000", \
     './dat/output_2D_B&W_MHD_y=x.dat' index (1) u 15:3 t "" w points pt 5 \
	ps 0.1

set ylabel '{/:Italic v_{/Symbol \136}}' font ",16"
set yrange [-0.3:0.7]
plot './dat/output_2D_B&W_MHD_y=x.dat' index (0) u 15:13 t "" w lines \
	 lc rgb "#000000", \
     './dat/output_2D_B&W_MHD_y=x.dat' index (1) u 15:13 t "" w points pt 5 \
	ps 0.1

set ylabel '{/:Italic v_{/Symbol \174\174}}' font ",16"
set yrange [-1.8:0.1]
plot './dat/output_2D_B&W_MHD_y=x.dat' index (0) u 15:14 t "" w lines \
	lc rgb "#000000", \
     './dat/output_2D_B&W_MHD_y=x.dat' index (1) u 15:14 t "" w points pt 5 \
	ps 0.1

set ylabel 'p_T' font ",16"
set yrange [0.6:1.9]
plot './dat/output_2D_B&W_MHD_y=x.dat' index (0) u 15:7 t "" w lines \
	lc rgb "#000000", \
     './dat/output_2D_B&W_MHD_y=x.dat' index (1) u 15:7 t "" w points pt 5 \
	ps 0.1

set xlabel '{/:Italic ||x+y||}' font ",16"
set ylabel '{/:Italic B_{/Symbol \136}}' font ",16"
set yrange [-1.2:1.2]
plot './dat/output_2D_B&W_MHD_y=x.dat' index (0) u 15:16 t "" w lines \
	lc rgb "#000000", \
     './dat/output_2D_B&W_MHD_y=x.dat' index (1) u 15:16 t "" w points pt 5 \
	ps 0.1

unset multiplot
unset xlabel

#--------------------------------------------------------------------------------
# Print the 5th Test Sod with cylindrical explosion

set terminal postscript eps enhanced color

set size ratio 1/1

unset xrange
unset yrange

set ticslevel 0.5

set pm3d at bs
set style fill transparent solid 0.8 noborder
set hidden3d
set palette

set parametric
set isosamples 51, 51

#set multiplot layout 4,1 rowsfirst

set output 'results/2D_Cyl_Expl_rho.eps'

set xlabel '{/:Italic x}' font ",16"
set ylabel '{/:Italic y}' font ",16"
set zlabel '{/:Italic {/Symbol r}}' font ",16"
set zrange [0.0:1.0]
splot './dat/output_2D_Cyl_Expl.dat' u 1:2:3 w pm3d notitle

set output 'results/2D_Cyl_Expl_v.eps'

set xlabel '{/:Italic x}' font ",16"
set ylabel '{/:Italic y}' font ",16"
set zlabel '||{/:Bold v}||' font ",16"
set zrange [-0.1:2.0]
splot './dat/output_2D_Cyl_Expl.dat' u 1:2:4 w pm3d notitle

set output 'results/2D_Cyl_Expl_p.eps'

set xlabel '{/:Italic x}' font ",16"
set ylabel '{/:Italic y}' font ",16"
set zlabel 'p' font ",16"
set zrange [0.0:1.1]
splot './dat/output_2D_Cyl_Expl.dat' u 1:2:7 w pm3d notitle

set output 'results/2D_Cyl_Expl_e.eps'

set xlabel '{/:Italic x}' font ",16"
set ylabel '{/:Italic y}' font ",16"
set zlabel '{/:Italic {/Symbol e}}' font ",16"
set zrange [1.7:2.95]
splot './dat/output_2D_Cyl_Expl.dat' u 1:2:11 w pm3d notitle

#unset multiplot
unset xlabel

#--------------------------------------------------------------------------------
# Print the 9th Test Orzag-

set terminal postscript eps enhanced color

set size ratio 1/1

unset xrange
unset yrange
unset zrange

set ticslevel 0.5

set view map

set pm3d at b map
set style fill transparent solid 0.8 noborder
set hidden3d

set palette gray

set parametric
set isosamples 51, 51

#set multiplot layout 4,1 rowsfirst

set output 'results/2D_OT_Vort_rho.eps'

set xlabel '{/:Italic x}' font ",16"
set ylabel '{/:Italic y}' font ",16"
set zlabel '{/:Italic {/Symbol r}}' font ",16"
#set zrange [0.0:1.0]
splot './dat/output_2D_OT_Vort.dat' i (1) u 1:2:3 w pm3d notitle

#set output 'results/2D_Cyl_Expl_v.eps'
#
#set xlabel '{/:Italic x}' font ",16"
#set ylabel '{/:Italic y}' font ",16"
#set zlabel '||{/:Bold v}||' font ",16"
#set zrange [-0.1:2.0]
#splot './dat/output_2D_Cyl_Expl.dat' u 1:2:4 w pm3d notitle
#
#set output 'results/2D_Cyl_Expl_p.eps'
#
#set xlabel '{/:Italic x}' font ",16"
#set ylabel '{/:Italic y}' font ",16"
#set zlabel 'p' font ",16"
#set zrange [0.0:1.1]
#splot './dat/output_2D_Cyl_Expl.dat' u 1:2:7 w pm3d notitle
#
#set output 'results/2D_Cyl_Expl_e.eps'
#
#set xlabel '{/:Italic x}' font ",16"
#set ylabel '{/:Italic y}' font ",16"
#set zlabel '{/:Italic {/Symbol e}}' font ",16"
#set zrange [1.7:2.95]
#splot './dat/output_2D_Cyl_Expl.dat' u 1:2:11 w pm3d notitle

#unset multiplot
unset xlabel


