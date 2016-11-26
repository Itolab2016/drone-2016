set grid
plot "flightlog.txt" u 1:6 w l title "ax" , "flightlog.txt" u 1:7 w l title "ay" , "flightlog.txt" u 1:8 w l title "az"
set terminal postscript eps  enhanced color
set output "accelerometer.eps"
replot
