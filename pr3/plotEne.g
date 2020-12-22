set terminal png size 800,600
set output "energy.png"
set title "Energy over time"
set xlabel "time [s]"
set ylabel "Energy"
plot "thermodynamics.dat" using 1:2 with lines t "Potential", \
"" using 1:3 with lines t "Kinetic", \
"" using 1:4 with lines t "Total"
pause -1