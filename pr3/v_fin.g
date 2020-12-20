set terminal png size 800,600
set output "v_fin.png"
binwidth=1
stats 'vel.dat'
bin(x,width)=width*floor(x/width)
set xrange [0:40]
plot 'vel.dat' using (bin($2,binwidth)):(1./STATS_records) smooth freq with boxes title "Final velocity distribution"
unset label
pause -1