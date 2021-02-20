set terminal png size 800,600
set output "2a_variances.png"
set title "Variances"
set xlabel "Longitut del bloc"
set ylabel "Variança"
set key top left
plot "thermodynamic_variances.dat" using 1:2 w l t "Variança de l'energia potencial", \
"" using 1:3 w l t "Variança de l'energia cinètica"

set key top right
set terminal png size 800,600
set output "2a_momentum.png"
set title "Total momentum over time"
set xlabel "temps [s]"
set ylabel "Moment"
plot "thermodynamics.dat" using 1:6 with lines title "Moment total"

unset yrange
set terminal png size 800,600
set output "2a_temperature.png"
set title "Evolució temporal de la temperatura"
set xlabel "temps [s]"
set ylabel "Temperatura [k]"
set key bottom left
plot "thermodynamics.dat" using 1:5 with lines title "Temperatura"

set key top right
set terminal png size 800,600
set output "2a_g(r)_ini.png"
set title "Distribució radial"
set xlabel "r [Å]"
set ylabel "Probabilitat"
plot "g(r)_ini.dat" using 1:2 with lines title "Configuració SC"

set terminal png size 800,600
set output "2a_g(r)_isolated.png"
set title "Distribució radial"
set xlabel "r [Å]"
set ylabel "Probabilitat"
plot "g(r)_isolated.dat" using 1:2 with lines title "Configuració inicial"

set terminal png size 800,600
set output "2a_g(r)_melted.png"
set title "Distribució radial"
set xlabel "r [Å]"
set ylabel "Probabilitat"
plot "g(r)_melted.dat" using 1:2 with lines title "Configuració final"

set terminal png size 800,600
set output "2a_ene.png"
set title "Evolució temporal de la energia"
set xlabel "temps [s]"
set ylabel "Energia"
plot "thermodynamics.dat" using 1:2 with lines title "Energia potencial", "" using 1:3 with lines title "Energia cinètica", "" using 1:4 with lines title "Energia total"
