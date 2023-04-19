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
unset ztics
set xrange [-0.5:0.5]
set yrange [-0.5:0.5]
set zrange [-0.5:0.5]
set output 'snapshot_random.png'
splot 'snapshot_random' with circles lw 2 lc rgb "black"
set output 'snapshot_minpos.png'
splot 'snapshot_minpos' with circles lw 2 lc rgb "black"
set output 'snapshot_mcequi.png'
splot 'snapshot_mcequi' with circles lw 2 lc rgb "black"
set output 'snapshot_mdinit.png'
splot 'snapshot_mdinit' with circles lw 2 lc rgb "black"
set output 'snapshot_mdequi.png'
splot 'snapshot_mdequi' with circles lw 2 lc rgb "black"
do for [i=1:5] {
  set output sprintf('snapshot_mc%04d.png', i)
  splot sprintf('snapshot_mc%04d', i) with circles lw 2 lc rgb "black"
  set output sprintf('snapshot_md%04d.png', i)
  splot sprintf('snapshot_md%04d', i) with circles lw 2 lc rgb "black"
}

set size nosquare
