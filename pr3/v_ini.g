set terminal png size 800,600
set output "v_ini.png"
binwidth=.1
stats 'vel.dat'
bin(x,width)=width*floor(x/width)
set xrange [0:2]
plot 'vel.dat' using (bin($1,binwidth)):(1./STATS_records) smooth freq with boxes title "Initial velocity distribution"
unset label
pause -1