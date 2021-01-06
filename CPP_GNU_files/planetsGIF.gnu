#! /usr/local/bin/gnuplot
reset
set terminal gif large animate size 900,900 font 'Verdana,10'
set output "planetsmotion_cartesian.gif"
n=150    #n frames
set encoding utf8

set view equal xy
set size square 
set tmargin 10
set rmargin 10
set bmargin 5
set lmargin 16

# Axes
set style line 11 lc rgb '#424141' lt 1
set border 3 back ls 11
set tics nomirror out scale 1
set xtics (-9.56245, -19.2195, -30.0798, 0, 9.56245, 19.2195, 30.0798)
set ytics (-5.21225, -9.56245, -19.2195, -30.0798, 0, 5.21225, 9.56245, 19.2195, 30.0798)
set xlabel 'Distance in A.U.'
set ylabel 'Distance in A.U.'
set arrow from graph 1,0 to graph 1.05,0 size screen 0.025,15,60 filled ls 11
set arrow from graph 0,1 to graph 0,1.05 size screen 0.025,15,60 filled ls 11
set xlabel 'Distance in A.U.'
set ylabel 'Distance in A.U.'

# Grid
set style line 12 lc rgb'#424141' lt 0 lw 1
set grid back ls 12

# Color definitions
set border linewidth 1.5
set style line 1 lc rgb '#FF6600' lt 7 lw 2		# Sun
set style line 2 lc rgb '#4A99A1' lt 7 lw 2		# Mercury
set style line 3 lc rgb '#996600' lt 7 lw 2 	# Venus
set style line 4 lc rgb '#000066' lt 7 lw 2		# Earth
set style line 5 lc rgb '#990000' lt 7 lw 2 	# Mars
set style line 6 lc rgb '#336600' lt 7 lw 2		# Jupiter
set style line 7 lc rgb '#5CC9F5' lt 7 lw 2 	# Saturn
set style line 8 lc rgb '#6638F0' lt 7 lw 2		# Uranus
set style line 9 lc rgb '#F4169F' lt 7 lw 2		# Neptun

set view equal xy
set size square 
set xrange [-32:32]
set yrange [-32:32]
set key at 35,35

set object circle at 0,0 size 5.21225 back fc '#C1C1C2' fs
set object circle at 0,0 size 9.56245 back fc '#C1C1C2' fs
set object circle at 0,0 size 19.2195 back fc '#C1C1C2' fs
set object circle at 0,0 size 30.0798 back fc '#C1C1C2' fs

i=0
load "animatePlanets.gnuplot"
set output
