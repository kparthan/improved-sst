# Gnuplot script file for plotting data in file "data"

set terminal post eps
set autoscale	# scale axes automatically
set xtic auto	# set xtics automatically
set ytic auto	# set ytics automatically
#set title ""
set xlabel "# of components"
set ylabel "message length (in bits)"
set output "msglens_components.eps"
plot "msglens_components.dat" using 1:5 notitle with points lc rgb "red"
