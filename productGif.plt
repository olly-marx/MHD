set terminal gif animate delay 10
set output filename.'.gif'
stats './dat/'.filename.'.dat' nooutput

do for [i=1:int(STATS_blocks)] {
    plot './dat/'.filename.'.dat' index (i-1) with lines
}
