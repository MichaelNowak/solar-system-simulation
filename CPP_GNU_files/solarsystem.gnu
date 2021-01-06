#! /usr/local/bin/gnuplot
reset
set terminal pngcairo size 900,900 enhanced font 'Verdana,10'
set encoding utf8
set polar
set grid polar
set angles radians
#unset key
set view equal xy
set size square 
set notics
set noytics
#set rtics (0.00465, 0.0093, 0.3871, 0.7233, 1, 1.5237, 5.2029, 9.5367, 19.1891, 30.0698) format ''
set rtics (0.00465, 0.0093, 0.396095, 0.724452, 1.00169, 1.53271, 5.21225, 9.56245, 19.2195, 30.0798) format ''
set tmargin 5
set rmargin 10
set bmargin 5
set lmargin 16

# Color definitions
set border linewidth 1.5
set style line 1 lc rgb '#931DD0' lt 1 lw 0.05
set style line 2 lc rgb '#4A99A1' lt 2 lw 1
set style line 3 lc rgb '#996600' lt 2 lw 1 
set style line 4 lc rgb '#000066' lt 2 lw 1
set style line 5 lc rgb '#990000' lt 2 lw 1 
set style line 6 lc rgb '#336600' lt 2 lw 1
set style line 7 lc rgb '#5CC9F5' lt 2 lw 1 
set style line 8 lc rgb '#6638F0' lt 2 lw 1
set style line 9 lc rgb '#990066' lt 2 lw 1

# Axes
set style line 11 lc rgb '#424141' lt 1
set border 3 back ls 11
set tics nomirror out scale 1
set xlabel 'Distance in A.U.'
set ylabel 'Distance in A.U.'
set arrow from graph 1,0 to graph 1.05,0 size screen 0.025,15,60 filled ls 11
set arrow from graph 0,1 to graph 0,1.05 size screen 0.025,15,60 filled ls 11
set xlabel 'Distance in A.U.'
set ylabel 'Distance in A.U.'

# Grid
set style line 12 lc rgb'#424141' lt 0 lw 1
set grid back ls 12

set xtics (0, 0.00465, 0.0093)
set ytics (0, 0.00465, 0.0093)
set rrange [0:0.01]
set output "SunMovement.png"
set object circle at 0,0 size 0.00465 fc '#F1E22A' fs transparent solid 0.2 noborder behind
plot 'SunOutput.txt' using 2:3 with linespoint title 'Sun' ls 1


set xtics (0, 0.396095, 1.00169, 1.53271)
set ytics (0, 0.396095, 0.724452, 1.00169, 1.53271)
set rrange [0:1.6]
set output "MeVEMaMovement.png"
plot 'MercuryOutput.txt' using 2:3 title 'Mercury' ls 2, 'VenusOutput.txt' using 2:3 title 'Venus' ls 3, 'EarthOutput.txt' using 2:3 title 'Earth' ls 4, 'MarsOutput.txt' using 2:3 title 'Mars' ls 5


set xtics (0, 9.56245, 19.2195, 30.0798)
set ytics (0, 5.21225, 9.56245, 19.2195, 30.0798)
set rrange [0:30.1]
set output "JSUNMovement.png"
plot 'JupiterOutput.txt' using 2:3 title 'Jupiter' ls 6, 'SaturnOutput.txt' using 2:3 title 'Saturn' ls 7, 'UranusOutput.txt' using 2:3 title 'Uranus' ls 8, 'NeptunOutput.txt' using 2:3 title 'Neptun' ls 9


reset
set terminal pngcairo size 900,900 enhanced font 'Verdana,10'
set output "partialMercuryMovement.png"
set view equal xy
set size square
set tmargin 5
set rmargin 10
set bmargin 5
set lmargin 13

# Color definitions
set border linewidth 1.5
set style line 1 lc rgb '#FF68D0' lt 1 lw 2 # --- pink
set style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 7 # --- red
set style line 3 lc rgb '#4A99A1' lt 2 lw 1

# Axes
set style line 11 lc rgb '#424141' lt 1
set border 3 back ls 11
set tics nomirror out scale 1
set xlabel 'Distance in A.U.'
set ylabel 'Distance in A.U.'
set arrow from graph 1,0 to graph 1.05,0 size screen 0.025,15,60 filled ls 11
set arrow from graph 0,1 to graph 0,1.05 size screen 0.025,15,60 filled ls 11
set xlabel 'Distance in A.U.'
set ylabel 'Distance in A.U.'

# Grid
set style line 12 lc rgb'#424141' lt 0 lw 1
set grid back ls 12

set trange [0:2*pi]

set parametric

# Radius
r = 0.396095
h = r / sqrt(2.)

set xrange [h-0.03:h+0.03]
set yrange [h-0.03:h+0.03]

set arrow from graph 0,0 to h,h ls 2 #nohead ls 2
set label 'r = 0.396095' at 0.265,0.26 textcolor ls 2

# Parametric functions for a circle
fx(t) = r*cos(t)
fy(t) = r*sin(t)

plot 'MercuryOutput.txt' using 4:5 with linespoint title 'Mercury' ls 3, fx(t),fy(t) ls 1 title 'Circle'

reset
set terminal pngcairo size 900,900 enhanced font 'Verdana,10'

set tmargin 5
set rmargin 10
set bmargin 5
set lmargin 16

# Grid
set style line 12 lc rgb'#424141' lt 0 lw 1
set grid back ls 12

# Axes
set style line 11 lc rgb '#424141' lt 1
set border 3 back ls 11
set border linewidth 1.5
set tics nomirror out scale 1

set xlabel 't [years]'

set ylabel 'Distance in A.U.'
set output "radiusFluctuationsMercuryZoom.png"
plot 'MercuryOutput.txt' every ::1000::1750 using 1:8 with linespoint title 'radius Mercury' ls 3

set output "radiusFluctuationsMercury.png"
plot 'MercuryOutput.txt' using 1:8 with linespoint title 'radius Mercury' ls 3

set ylabel 'Velocity in A.U./yr'
set output "velocityFluctuationsMercuryZoom.png"
plot 'MercuryOutput.txt' every ::1000::1750 using 1:9 with linespoint title 'velocity Mercury' ls 3

set output "velocityFluctuationsMercury.png"
plot 'MercuryOutput.txt' using 1:9 with linespoint title 'velocity Mercury' ls 3