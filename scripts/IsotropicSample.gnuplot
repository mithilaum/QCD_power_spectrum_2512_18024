#!/usr/bin/gnuplot --persist

# Set up PDF output using cairolatex terminal with specific styling
set term cairolatex pdf size 3 in, 2 in \
	standalone \
    dashed \
    linewidth 1.5 \
    mono \
    rounded \
    header "\\usepackage{amsmath}\n"
   
# Handle missing data points
set datafile missing "0.000000e+00"

# Function to handle non-negative values (corrects small numerical errors)
# Uses a very small positive number (1e-9) for non-positive values
makeNN(x) = ((x <= 0.) ? 1e-9 : x)

# Set margins for plot layout
set lmargin at screen 0.21
set rmargin at screen 0.95
set tmargin at screen 0.88
set bmargin at screen 0.20

# Define line styles
set linestyle 11 lt 1 lw 1 lc 0
set linestyle 12 lt 3 lw 1 lc 9
set linestyle 13 pt 19 ps 0.7 lw 3 lc 0 

# Set axis labels with LaTeX formatting
set xlabel '$\ell$' offset 0, 0.2
set ylabel '$H_\ell^{}$' offset -0.5

# Configure y-axis scaling and style
set logscale y
set style data lines
set xrange [:70]
set format y '$10^{%L}$'

# Set up axis ticks
set ytics 1e-8, 100, 1 add ('1' 1) nomirror
set mxtics 10
set mytics 10
set xtics nomirror
set x2tics -10, 2, -8 add ('$\;\,3^\circ$' 120, '$\;\,6^\circ$' 60, \
	'$\;\,12^\circ$' 30, '$\;\,24^\circ$' 15, '$\;\;180^\circ$' 2)

# Define asymptotic value
H_asym = 1. + 1./3.592

# First plot: Tracks
set output 'Hl_isotropic_tracks.tex'
set yrange [1e-9:1.]

# Add labels for different N* values
# set label 1 at 24.5, 2.8e-2 '{\footnotesize $N^*=15$}'
set label 2 at 10, 3.25e-2 '{\footnotesize $N^*=126$}'
set label 3 at 27, 2.7e-3 '{\footnotesize $N^*=1{,}020$}'
set label 4 at 42, 8.2e-6 '{\footnotesize $N^*=8{,}146$}'
set label 5 at 52, 1e-8 '{\footnotesize $N^*=19{,}891$}'

# Configure secondary y-axis
set logscale y2
set y2range [1e-9:1.]
set y2tics 10, 10, 20 out add ('' 9.695e-03, '' 1.276e-03, '' 1.558e-04, '' 5.0274e-05)


# Plot tracks data
plot   	"./dat/powerspectrum_detector_holes_etaMax_3.009_normalized.dat" using 1:2 with points pt 6 ps 0.35 lc rgb "#FA8072"  notitle, \
	"./dat/128_6.dat" using 1:(makeNN($2)) notitle ls 11, \
	H_asym / 126 notitle ls 12, \
	"./dat/1024_6.dat" using 1:(makeNN($2)) notitle ls 11, \
	H_asym / 1020 notitle ls 12, \
	"./dat/8192_6.dat" using 1:(makeNN($2)) notitle ls 11, \
	H_asym / 8146 notitle ls 12, \
	"./dat/H_ISO_thomson_20000_6_3.009000_3.009000.dat" using 1:(makeNN($2)) notitle ls 11, \
	H_asym / 19891 notitle ls 12


# Second plot: Towers
set output 'Hl_isotropic_towers.tex'

# Configure x-axis for towers plot
unset logscale x
unset logscale x2
set xrange [:70]

# Configure y-axis for towers plot
set yrange [1e-7:1.]

# Add vertical arrows
set arrow 1 from 5.8, 1e-7  to 5.8, 8.5e-2 nohead ls 12
set arrow 2 from 24.6, 1e-7 to 24.6, 9e-3  nohead ls 12 
set arrow 3 from 58.4, 1e-7 to 58.4, 2e-3  nohead ls 12
set arrow 4 from 60.3, 1e-7 to 60.3, 2e-3  nohead ls 12
# CM: Arrows 2 and 4 originally had lc 9. 

# Configure y2 axis for towers plot
set y2range [1e-7:1.]
unset y2tics
set y2tics 10, 10, 20 out add ('' 1.035e-01, '' 9.695e-03, '' 1.467e-03, '' 8.991e-04) 

# Update labels for towers plot
set label 1 at 23, 3e-2 '{\footnotesize $N^*=15$}'
set label 2 at 29, 3e-3 '{\footnotesize $N^*=126$}'
set label 3 at 35, 2e-4 '{\footnotesize $N^*=818$}'
set label 4 at 40, 7e-6 '{\footnotesize $N^*=1{,}144$}'
unset label 5

# Plot towers data
plot "./dat/16_6.dat" using 1:4 notitle ls 11, \
	"./dat/128_6.dat" using 1:4 notitle ls 11, \
	"./dat/1024_6.dat" using 1:4 notitle ls 11, \
	"./dat/8192_6.dat" using 1:4 notitle ls 11
	