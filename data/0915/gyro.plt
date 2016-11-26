set grid
plot "flightlog.txt" u 1:3 w l title "gx" ,"flightlog.txt" u 1:4 w l title "gy" , "flightlog.txt" u 1:5 w l title "gz"
set terminal postscript eps  enhanced color
set output "gyro.eps"
replot
