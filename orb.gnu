#!/usr/bin/gnuplot
set term gif enhanced animate

set output 'Orb_5.gif'
set title "{/Symbol Y}(x,t)"

set xrange [-100:100]
#set yrange [-1.6:1.6]
#set xrange [0:14000e-10]
#set yrange [0:14000e-10]
#set zrange [0:14000e-10]
set xlabel "Position (Bohr Radii)"
set ylabel "Amplitude"

i=0
n=999

load "orb.gnuplot"
