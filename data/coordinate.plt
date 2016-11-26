set grid
plot "flightlog.txt" u 1:26 w l title "x" , "flightlog.txt" u 1:27 w l title "y" , "flightlog.txt" u 1:28 w l title "z"
set terminal postscript eps  enhanced color
set output "coordinate.eps"
replot
