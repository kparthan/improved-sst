# Gnuplot script file for plotting data in file "data"

set terminal post eps
set autoscale	# scale axes automatically
set xtic auto	# set xtics automatically
set ytic auto	# set ytics automatically
set xr [0:]
set xlabel "# of components"
set ylabel "message length (in bits)"
set output "/home/pkas7/Research/Work/improved-sst/mixture/simulation/msglens-infer.eps"
plot "/home/pkas7/Research/Work/improved-sst/mixture/simulation/msglens-infer.dat" using 1:2 notitle with linespoints lc rgb "red"
