# 1:now_time 2:interval
# 3:gx 4:gy 5:gz 6:ax 7:ay 8:az 9:mx 10:my 11:mz
# 12:roll 13:pitch
# 14:R 15:L 16:F 17:B

set grid
plot "flightlog.txt" u 1:3 w l t "gx" ,"flightlog.txt" u 1:4 w l t "gy" , "flightlog.txt" u 1:5 w l t "gz"
set terminal postscript eps  enhanced color
set output "gyro.eps"
replot
