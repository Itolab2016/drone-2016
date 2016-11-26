set grid
plot "flightlog.txt" u 1:9 w l title "mx" , "flightlog.txt" u 1:10 w l title "my" , "flightlog.txt" u 1:11 w l title "mz"
set terminal postscript eps  enhanced color
set output "geomagnetic.eps"
replot
