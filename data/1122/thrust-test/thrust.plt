set grid
plot "thrust.txt" u 1:2 w lp lw 5 pt 6 title "thrust"
set yrange [-20:300]
set xlabel 'pwm[ms]'
set ylabel 'thrust[g]'
set terminal postscript eps  enhanced color
set output "thrust.eps"
replot
