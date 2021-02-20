set terminal png size 800,600
set output "energies.png"
set title "Energy vs density"
set xlabel "density [g/cm³]"
set ylabel "energy [kJ/mol]"
set key left top
plot "thermodynamics.dat" using 1:2 pt 7 ps 2 t "Potential", \
"" u 1:3 pt 7 ps 2 t "Kinetic", \
"" u 1:4 pt 7 ps 2 t "Total"

set terminal png size 800,600
set output "pressure.png"
set title "Energy vs density"
set xlabel "density [g/cm³]"
set ylabel "Pressure [Pa]"
plot "thermodynamics.dat" using 1:5 pt 7 ps 2 t "Pressure"

set terminal png size 800,600
set output "g(r).1.png"
set title "Funció de distribució per a ρ'=0.1"
set xlabel "Distncia [Å]"
set ylabel "Probabilitat"
plot "g(r).1.dat" using 1:2 w l t "Distribució radial"

set terminal png size 800,600
set output "g(r).2.png"
set title "Funció de distribució per a ρ'=0.2"
set xlabel "Distncia [Å]"
set ylabel "Probabilitat"
plot "g(r).2.dat" using 1:2 w l t "Distribució radial"

set terminal png size 800,600
set output "g(r).4.png"
set title "Funció de distribució per a ρ'=0.1"
set xlabel "Distncia [Å]"
set ylabel "Probabilitat"
plot "g(r).4.dat" using 1:2 w l t "Distribució radial"

set terminal png size 800,600
set output "g(r).6.png"
set title "Funció de distribució per a ρ'=0.6"
set xlabel "Distncia [Å]"
set ylabel "Probabilitat"
plot "g(r).6.dat" using 1:2 w l t "Distribució radial"

set terminal png size 800,600
set output "g(r).8.png"
set title "Funció de distribució per a ρ'=08"
set xlabel "Distncia [Å]"
set ylabel "Probabilitat"
plot "g(r).8.dat" using 1:2 w l t "Distribució radial"