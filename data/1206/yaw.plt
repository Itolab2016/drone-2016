# 1:now_time 2:interval
# 3:gx 4:gy 5:gz 6:ax 7:ay 8:az 9:mx 10:my 11:mz
# 12:roll 13:pitch 14:yaw
# 15:R 16:L 17:F 18:B

set grid
plot "flightlog.txt" u 1:14 w l t "yaw"
set terminal postscript eps  enhanced color
set output "yaw.eps"
replot
