set grid;
plot "flightlog.txt" u 10:1 with lines lw 5 title "U";
set terminal postscript eps enhanced color;
set output "U.eps";
replot;

set grid;
plot "flightlog.txt" u 10:2 with lines lw 5 title "V";
set terminal postscript eps enhanced color;
set output "V.eps";
replot;

set grid;
plot "flightlog.txt" u 10:3 with lines lw 5 title "W";
set terminal postscript eps enhanced color;
set output "W.eps";
replot;

set grid;
plot "flightlog.txt" u 10:4 with lines lw 5 title "p";
set terminal postscript eps enhanced color;
set output "p.eps";
replot;

set grid;
plot "flightlog.txt" u 10:5 with lines lw 5 title "q";
set terminal postscript eps enhanced color;
set output "q.eps";
replot;

set grid;
plot "flightlog.txt" u 10:6 with lines lw 5 title "r";
set terminal postscript eps enhanced color;
set output "r.eps";
replot;

set grid;
plot "flightlog.txt" u 10:7 with lines lw 5 title "phi";
set terminal postscript eps enhanced color;
set output "phi.eps";
replot;

set grid;
plot "flightlog.txt" u 10:8 with lines lw 5 title "th";
set terminal postscript eps enhanced color;
set output "th.eps";
replot;

set grid;
plot "flightlog.txt" u 10:9 with lines lw 5 title "yaw";
set terminal postscript eps enhanced color;
set output "yaw.eps";
replot;
