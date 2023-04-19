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

set terminal png enhanced
set size 1.0,1.0
set size square
unset xtics
unset ytics
set xrange [-0.5:0.5]
set yrange [-0.5:0.5]
set output 'snapshot_random.png'
plot 'snapshot_random' with circles lw 2 lc rgb "black"
set output 'snapshot_minpos.png'
plot 'snapshot_minpos' with circles lw 2 lc rgb "black"
set output 'snapshot_mcequi.png'
plot 'snapshot_mcequi' with circles lw 2 lc rgb "black"
set output 'snapshot_mdinit.png'
plot 'snapshot_mdinit' with circles lw 2 lc rgb "black"
set output 'snapshot_mdequi.png'
plot 'snapshot_mdequi' with circles lw 2 lc rgb "black"
do for [i=1:5] {
  set output sprintf('snapshot_mc%04d.png', i)
  plot sprintf('snapshot_mc%04d', i) with circles lw 2 lc rgb "black"
  set output sprintf('snapshot_md%04d.png', i)
  plot sprintf('snapshot_md%04d', i) with circles lw 2 lc rgb "black"
}

set size nosquare
