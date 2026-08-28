# Gnuplot plotting script used by plot_live_data script

#
# Purpose: Plot energies and growth rates during or after a JOREK code run
#
# Date: 2011-03-30
# Author: Matthias Hoelzl, IPP Garching
#

no0=<no0>
ps=<ps>
qtty="<qtty>"
title="<title>"
ncols0=2+no0
ncols=<ncols>
logy=<logy>
xlabel="<xlabel>"
ylabel="<ylabel>"
x_toSI=<qtty_x2si>
y_toSI=<qtty_y2si>

if ( ps==1 ) set term postscript enhanced color
if ( ps==1 ) set output qtty.'.ps'
if ( logy==1 ) set log y
set key outside
set title qtty.' '.title
set xlabel xlabel
set ylabel ylabel
set xrange <xrange>
set yrange <yrange>
set format y "%g"
set format x "%g"
set border lc 'white'
set ylabel tc 'white'
set xlabel tc 'white'
set key tc 'white'
set grid lc 'white'

plot for [i=ncols0:ncols+1] qtty.'.dat' u ($1*x_toSI):(column(i)*y_toSI) w lp lc i t columnhead(i)

if ( ps==0 ) print ''
if ( ps==0 ) print 'Hints:'
if ( ps==0 ) print '* Use right mouse button to zoom'
if ( ps==0 ) print '* Press "U" to unzoom again'
if ( ps==0 ) print '* Press "E" to update the data'
if ( ps==0 ) print '* Click with left mouse button to exit'
if ( ps==0 ) print ''
if ( ps==0 ) pause mouse button1
