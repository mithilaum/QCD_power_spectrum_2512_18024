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

# make non-negative (to correct small numerical errors in nearly zero moments)
makeNN(x) = ((x <= 0.) ? 0. : x)

# Pin the corners of the graph to the canvas
set lmargin at screen 0.18
set rmargin at screen 0.95

set tmargin at screen 0.95
set bmargin at screen 0.2

set xlabel '$\ell$' # offset -1.5, 0.7
set ylabel '$H_\ell^{}$' offset 1.

set linestyle 11 lt 1 lc 9 lw 1 pt 19 ps 0.4
set linestyle 12 lt 1 lc 0 lw 2 pt 7 ps 0.4
set linestyle 13 lt 1 lc 9 lw 1 pt 17 ps 0.4
set linestyle 14 lt 1 lc 0 lw 2 pt 18 ps 0.4
set linestyle 15 lt 2 lc 0 

# set style data linespoints
set style data lines
														
set output 'threeJet_3_asym.tex'
dat = './dat/H_ee_jjj_03_all_spherical_evenOdd.dat'

set xtics 0, 5, 200 #add ('' 5 0, '' 15 0, '' 25 0, '' 35 0, '' 45 0, '' 55 0, '' 65 0)
set mxtics 5

set key at screen 0.57, 0.9 spacing 1.2

set ytics nomirror
set y2tics -2, 1, -1 add ('' 7.48836e-02, '' 3.47579e-01) out scale 1.

set label 1 at 22.5, 0.48 '$n=3$'
set label 2 at 22, 0.12 '$N=60$'
				
plot [0:30][0.:1.05] 7.48836e-02 ls 15 title '$\langle f | f \rangle$', \
							3.47579e-01 ls 15 lc 0 notitle, \
							dat using 1:2 ls 11 title 'odd', \
							'' using ($1+1):(makeNN($3)) ls 12 title 'even', \
							'' using 1:($4) ls 11 notitle, \
							'' using ($1+1):(makeNN($5)) ls 12 notitle
# set style data linespoints
set style data lines

set output 'threeJet_78_asym.tex'
dat = './dat/H_ee_jjj_78_all_spherical_evenOdd.dat'

#~ set key at 80, 0.97 spacing 1.2

set xtics 0, 25, 200
set mxtics 5

unset y2tics
set y2tics -2, 1, -1 add ('' 1.58165e-01, '' 4.25262e-01) out scale 1.

set label 1 at 130, 0.56 '$n=3$'
set label 2 at 123, 0.22 '$N=57$'
														
plot [0:160][0.:1.05] 1.58165e-01 ls 15 title '$\langle f | f \rangle$', \
							4.25262e-01 ls 15 lc 0 notitle, \
							dat using 1:2 ls 11 title 'odd', \
							'' using ($1+1):(makeNN($3)) ls 12 title 'even', \
							'' using 1:($4) ls 11 notitle, \
							'' using ($1+1):(makeNN($5)) ls 12 notitle
