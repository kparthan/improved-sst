# Gnuplot script file for plotting data in file "data"

set terminal post eps
set autoscale	# scale axes automatically
set xtic auto	# set xtics automatically
set ytic auto	# set ytics automatically
set title "# of components: 250"
set xlabel "# of iterations"
set ylabel "message length (in bits)"
set output "/home/pkas7/Research/Work/improved-sst/mixture/plots/normal_weights_update/250.eps"
plot "/home/pkas7/Research/Work/improved-sst/mixture/msglens/normal_weights_update/250.dat" using 1:2 notitle with linespoints lc rgb "red"
