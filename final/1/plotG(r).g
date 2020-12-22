set terminal png size 800,600
set output "g(r).png"
set title "g(r) crystal"
set xlabel "distance [Å]"
set ylabel "Probability [%]"
plot "g(r)_melted.dat" using 1:2 w l t "After melting", \
"g(r)_isolated.dat" using 1:2 w l t "After isolating"

set terminal png size 800,600
set output "g(r)_ini.png"
set title "g(r) crystal"
set xlabel "distance [Å]"
set ylabel "Probability [%]"
plot "g(r)_ini.dat" using 1:2 w l t "Initial configuration"