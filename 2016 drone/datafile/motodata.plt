set grid;
plot "flightlog.txt" u 10:1 w l title "U";
set terminal postscript eps enhanced color;
set output "U.eps";
replot;

set grid;
plot "flightlog.txt" u 10:2 w l title "V";
set terminal postscript eps enhanced color;
set output "V.eps";
replot;

set grid;
plot "flightlog.txt" u 10:3 w l title "W";
set terminal postscript eps enhanced color;
set output "W.eps";
replot;

set grid;
plot "flightlog.txt" u 10:4 w l title "p";
set terminal postscript eps enhanced color;
set output "p.eps";
replot;

set grid;
plot "flightlog.txt" u 10:5 w l title "q";
set terminal postscript eps enhanced color;
set output "q.eps";
replot;

set grid;
plot "flightlog.txt" u 10:6 w l title "r";
set terminal postscript eps enhanced color;
set output "r.eps";
replot;

set grid;
plot "flightlog.txt" u 10:7 w l title "phi";
set terminal postscript eps enhanced color;
set output "phi.eps";
replot;

set grid;
plot "flightlog.txt" u 10:8 w l title "th";
set terminal postscript eps enhanced color;
set output "th.eps";
replot;

set grid;
plot "flightlog.txt" u 10:9 w l title "yaw";
set terminal postscript eps enhanced color;
set output "yaw.eps";
replot;
