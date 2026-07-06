#!/bin/bash

gnuplot -persist <<-EOFMarker
    set title 'phase'
    set xlabel 'x'
    set ylabel 'f'
	plot 'f_prof.dat' using 1:2 with lines title 'y=25%', '' using 1:3 with lines title 'y=50%', '' using 1:4 with lines title 'y=75%'
EOFMarker

gnuplot -persist <<-EOFMarker
    set title 'velocity in x'
    set xlabel 'u_x'
    set ylabel 'y'
	plot 'u_x_prof.dat' using 2:1 with lines title 'x=25%', '' using 3:1 with lines title 'x=50%', '' using 4:1 with lines title 'x=75%'
EOFMarker

gnuplot -persist <<-EOFMarker
    set title 'pressure'
    set xlabel 'x'
    set ylabel 'p'
	plot 'p_prof.dat' using 1:2 with lines title 'y=25%', '' using 1:3 with lines title 'y=50%', '' using 1:4 with lines title 'y=75%'
EOFMarker

