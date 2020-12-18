set title "Temperature over time"
set xlabel "time [s]"
set ylabel "Temperature"
plot "temp.dat" using 1:2 with lines title "Temperature over time"
pause -1