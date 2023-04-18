# Gnuplot script file for plotting data
set encoding utf8
set autoscale                                                       # scale axes automatically
unset logscale                                                      # remove any log-scaling
unset xlabel                                                        # remove any previous xlabel
unset ylabel                                                        # remove any previous ylabel
set xtics auto                                                      # set xtics automatically
set ytics auto                                                      # set ytics automatically
unset key
unset size

cd sprintf('%s\%.3f\%.3f', path, temp, dens)
set terminal postscript eps enhanced color
set size 0.5,0.5
set xlabel 'r^*'
set ylabel 'g(r^*)'
set xrange [0:6]
set output 'mc_corr.eps'
plot 'mc_corr' with lines lc rgb "black"
set output 'md_corr.eps'
plot 'md_corr' with lines lc rgb "black"

cd path