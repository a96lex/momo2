set terminal png size 800,600
set output "variances.png"
set title "Variances"
set xlabel "Block Length"
set ylabel "Variance"
plot "thermodynamic_variances.dat" using 1:2 w l t "Pot var", \
"" using 1:3 w l t "Kin var"
pause -1