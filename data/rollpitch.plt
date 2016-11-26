set grid
plot "flightlog.txt" u 1:12 w l title "roll" , "flightlog.txt" u 1:13 w l title "pitch"
set xlabel 'time[us]'
set ylabel 'value[deg]'
set terminal postscript eps  enhanced color
set output "rp.eps"
replot
