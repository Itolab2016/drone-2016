set grid
plot "flightlog.txt" u 1:22 w l title "R" , "flightlog.txt" u 1:23 w l title "L" , "flightlog.txt" u 1:24 w l title "F" , "flightlog.txt" u 1:25 w l title "B"
set yrange [0.5:2.5]
set xlabel 'time[us]'
set ylabel 'output[ms]'
set terminal postscript eps  enhanced color
set output "RLFB.eps"
replot
