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

set terminal gif enhanced animate delay 1 optimize
set size 1.0,0.5
unset xtics
unset ytics
halflbox=0.5
set xrange [-halflbox:-halflbox+0.2]
set yrange [-0.1:0.1]
set output 'movie.gif'
do for [i=1:100] {
  plot sprintf('frame_01%06d', i) u 1:(0):2 notitle with circles lw 2 lc rgb "black"
}

set size nosquare
