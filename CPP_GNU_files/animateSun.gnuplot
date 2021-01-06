# Plot every (i*1003) line (i*1 year)

set style circle radius 0.00465
plot 'SunOutput.txt' every ::0::(i*1003) using 4:5 w l lt 1 lw 1.5 title sprintf("t=%i",i), 'SunOutput.txt' every ::(i*1003)::(i*1003) using 4:5 with circles fc '#D62722' title 'Sun'
i=i+1
if (i < n) reread