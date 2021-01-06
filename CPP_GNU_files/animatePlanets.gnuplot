# Plot every i*1003 line (i*1 year)

plot 'MercuryOutput.txt' every ::(i*1003)::(i*1003) using 4:5 title sprintf("t=%i     Mercury",i) ls 2, 'VenusOutput.txt' every ::(i*1003)::(i*1003) using 4:5 title 'Venus' ls 3, 'EarthOutput.txt' every ::(i*1003)::(i*1003) using 4:5 title 'Earth' ls 4, 'MarsOutput.txt' every ::(i*1003)::(i*1003) using 4:5 title 'Mars' ls 5, 'JupiterOutput.txt' every ::(i*1003)::(i*1003) using 4:5 title 'Jupiter' ls 6, 'SaturnOutput.txt' every ::(i*1003)::(i*1003) using 4:5 title 'Saturn' ls 7, 'UranusOutput.txt' every ::(i*1003)::(i*1003) using 4:5 title 'Uranus' ls 8, 'NeptunOutput.txt' every ::(i*1003)::(i*1003) using 4:5 title 'Neptun' ls 9

i=i+1
if (i < n) reread