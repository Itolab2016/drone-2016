set grid;
plot "flightlog.txt" u 10:1 with lines lw 5 title "U","flightlog.txt" u 10:2 with lines lw 5 title "V","flightlog.txt" u 10:3 with lines lw 5 title "W","flightlog.txt" u 10:4 with lines lw 5 title "p","flightlog.txt" u 10:5 with lines lw 5 title "q","flightlog.txt" u 10:6 with lines lw 5 title "r","flightlog.txt" u 10:7 with lines lw 5 title "phi","flightlog.txt" u 10:8 with lines lw 5 title "th","flightlog.txt" u 10:9 with lines lw 5 title "yaw";
set terminal postscript eps enhanced color;
set output "alldata.eps";
replot;
