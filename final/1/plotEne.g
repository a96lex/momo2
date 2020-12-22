set title "Temperature over time"
set xlabel "time [s]"
set ylabel "Temperature"
plot "thermodynamics.dat" using 1:2 with lines title "Potential", "" using 1:3 with lines title "Kinetic", "" using 1:4 with lines title "Total"
pause -1