#!/usr/bin/gnuplot --persist

# Send the output to a tex file which can be compiled to pdf using ... [user]$ pdflatex dOmegaDx.tex
# Make the font sans-serif (including math), consistent with a beamer presentation
set term cairolatex pdf size 3 in, 2 in \
	standalone \
   dashed \
   linewidth 1.5 \
   mono \
   rounded \
   header "\\usepackage{amsmath}\n"
	#~ header "\\renewcommand{\\familydefault}{\\sfdefault}\n\\usepackage{sansmath}\n\\sansmath\n"

# Pin the corners of the graph to the canvas
set lmargin at screen 0.22
set rmargin at screen 0.95
set bmargin at screen 0.22
#~ set tmargin at screen 0.95

unset grid

set output 'Hl_attenuation.tex'

set style data lines

set xlabel '$\ell$'
set ylabel '$\bar{h}_\ell^{\,2}$' offset 0.

lMax = 80

set xtics 0, 10, lMax
set mxtics 5
set key at 79, 0.96

set ytics 0, 0.2, 1

plot [0:lMax][0:1] "./dat/hl.dat" \
		using 1:($2**2) ls 1 lw 1.5 title '$\lambda=1^\circ\phantom{.00}$', \
	'' using 1:($3**2) ls 2 lw 1.5 title '$\lambda=2.33^\circ$', \
	'' using 1:($4**2) ls 3 lw 1.5 title '$\lambda=6^\circ\phantom{.00}$'
	
set output 'Hl_attenuation_log.tex'

lMax = 460
set key at 460, 0.2

set ylabel offset 0.5
set lmargin at screen 0.23
set rmargin at screen 0.93

set logscale y
set ytics 1e-16, 1e4, 1 add ('1' 1)
set xtics 0, 100, 500
set format y '$10^{%L}$'

plot [0:lMax][1e-16:1] "./dat/hl.dat" \
		using 1:($2**2) ls 1 lw 1.5 title '$\lambda=1^\circ\phantom{.00}$', \
	'' using 1:($3**2) ls 2 lw 1.5 title '$\lambda=2.33^\circ$', \
	'' using 1:($4**2) ls 3 lw 1.5 title '$\lambda=6^\circ\phantom{.00}$'

