set grid
plot "flightlog.txt" u 1:14 w l title "yaw"
set xlabel 'time[us]'
set ylabel 'value[deg]'
set terminal postscript eps  enhanced color
set output "yaw.eps"
replot
