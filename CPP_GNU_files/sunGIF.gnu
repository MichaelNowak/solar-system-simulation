#! /usr/local/bin/gnuplot
reset
set terminal gif large animate size 900,900 font 'Verdana,10'
set output "sun.gif"
n=150    #n frames
set view equal xy
set size square 
set tmargin 5
set rmargin 10
set bmargin 5
set lmargin 16

# Axes
set style line 11 lc rgb '#424141' lt 1
set border 3 back ls 11
set tics nomirror out scale 1
set xtics (-0.0093, -0.00465, 0, 0.00465, 0.0093)
set ytics (-0.0093, -0.006975, -0.00465, -0.002325, 0, 0.002325, 0.00465, 0.006975, 0.0093)
set xlabel 'Distance in A.U.'
set ylabel 'Distance in A.U.'
set arrow from graph 1,0 to graph 1.05,0 size screen 0.025,15,60 filled ls 11
set arrow from graph 0,1 to graph 0,1.05 size screen 0.025,15,60 filled ls 11
set xlabel 'Distance in A.U.'
set ylabel 'Distance in A.U.'

# Grid
set style line 12 lc rgb'#424141' lt 0 lw 1
set grid back ls 12

set xrange [-0.015:0.015]
set yrange [-0.015:0.015]
set key at graph 1.03,1.03
set object circle at 0,0 size 0.00465 fc '#F1E22A' fs transparent solid 0.2 noborder behind


i=0
load "animateSun.gnuplot"
set output