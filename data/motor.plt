set grid
plot "motortest.txt" u 1:4 w l title "vibration" , "motortest.txt" u 1:5 w l title "output"
set terminal postscript eps  enhanced color
set output "motortest.eps"
replot
