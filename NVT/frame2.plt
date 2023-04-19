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
set size 1.0,1.0
set size square
unset xtics
unset ytics
set xrange [-0.5:0.5]
set yrange [-0.5:0.5]
set output 'movie.gif'
do for [i=1:100] {
  plot sprintf('frame_01%06d', i) notitle with circles lw 2 lc rgb "black"
}

set size nosquare
