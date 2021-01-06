#!/usr/bin/gnuplot

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

# Plot

set ylabel 'Momentum [(Earth masses) A.U./yr]'
set output 'momentum_cons.png'
plot "Energy_cons.txt"  using 1:5 title "p_{tot}"

set ylabel 'Angular momentum [(Earth masses) A.U.^{2}/yr]'
set output 'angular_momentum_cons.png'
plot "Energy_cons.txt"  using 1:6 title "L_{tot}"


set ylabel 'Energy [(Earth masses) A.U.^{2}/yr^{2}]'
set output 'energy_cons_kin.png'
plot "Energy_cons.txt"  using 1:2 title "E_{kin}"

set output 'energy_cons_pot.png'
plot "Energy_cons.txt"  using 1:3 title "E_{pot}"

set output 'energy_cons_tot.png'
plot "Energy_cons.txt"  using 1:4 title "E_{tot}"

set yrange [-2000:3500]
set output 'energy_cons.png'
plot "Energy_cons.txt"  using 1:2 title "E_{kin}", "Energy_cons.txt"  using 1:3 title "E_{pot}", "Energy_cons.txt"  using 1:4 title "E_{tot}"