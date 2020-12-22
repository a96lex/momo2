set terminal png size 800,600
set output "temperature.png"
set title "Temperature over time"
set xlabel "time [s]"
set ylabel "Temperature"
plot "thermodynamics.dat" using 1:5 with lines t "Temperature"
pause -1