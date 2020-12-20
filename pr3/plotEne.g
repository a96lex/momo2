set title "Temperature over time"
set xlabel "time [s]"
set ylabel "Temperature"
plot "thermodynamics.dat" using 1:5 with lines t "Temperature"
pause -1