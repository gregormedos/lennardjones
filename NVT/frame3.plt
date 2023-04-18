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
unset ztics
halflbox=0.5
set xrange [-halflbox:halflbox]
set yrange [-halflbox:halflbox]
set zrange [-halflbox:halflbox]
set output 'movie.gif'
do for [i=1:100] {
  splot sprintf('frame_01%06d', i) notitle with circles lw 2 lc rgb "black"
}

set size nosquare
