set title "Temperature over time"
set xlabel "time [s]"
set ylabel "Temperature"
plot "thermodynamicsF.dat" using 1:2 with lines title "Potential", "" using 1:3 with lines title "Kinetic", "" using 1:4 with lines title "Total", "" using 1:5 with lines title "Temperatura"
pause -1