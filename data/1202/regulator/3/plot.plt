# 1:now_time 2:interval
# 3:gx 4:gy 5:gz 6:ax 7:ay 8:az 9:mx 10:my 11:mz
# 12:roll 13:pitch 14:yaw 15:rroll 16:rpitch 17:ryaw
# 18:stickRx 19:stickRy 20:stickLx 21:stickLy
# 22:R 23:L 24:F 25:B
# 26:x 27:y 28z

set grid
plot "flightlog.txt" u 1:3 w l ,"flightlog.txt" u 1:4 w l , "flightlog.txt" u 1:5 w l , "flightlog.txt" u 1:6 w l , "flightlog.txt" u 1:7 w l , "flightlog.txt" u 1:8 w l , "flightlog.txt" u 1:9 w l , "flightlog.txt" u 1:10 w l , "flightlog.txt" u 1:11 w l , "flightlog.txt" u 1:12 w l , "flightlog.txt" u 1:13 w l , "flightlog.txt" u 1:14 w l , "flightlog.txt" u 1:15 w l , "flightlog.txt" u 1:16 w l , "flightlog.txt" u 1:17 w l , "flightlog.txt" u 1:18 w l , "flightlog.txt" u 1:19 w l , "flightlog.txt" u 1:20 w l , "flightlog.txt" u 1:21 w l , "flightlog.txt" u 1:22 w l , "flightlog.txt" u 1:23 w l , "flightlog.txt" u 1:24 w l , "flightlog.txt" u 1:25 w l , "flightlog.txt" u 1:26 w l , "flightlog.txt" u 1:27 w l , "flightlog.txt" w 1:28 w l
set terminal postscript eps  enhanced color
set output "data.eps"
replot
