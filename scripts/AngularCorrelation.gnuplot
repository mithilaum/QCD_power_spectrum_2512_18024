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
set lmargin at screen 0.13
set rmargin at screen 0.96

set bmargin at screen 0.22
set tmargin at screen 0.95

set style data lines

deg(x)=x*180./pi

set xlabel '$\cos\xi$' offset 0.5
set ylabel '$A(\cos\xi)$' offset 0.25

set key at 0.9, 8.5 spacing 1.1 opaque width -9 box
# set key at 0.62, 6.5 spacing 1.1 opaque width -9 box

yMin = 0.
yMax = 9.5
zMM = 1.03

set yrange [yMin:yMax]

set xtics -1., 0.5, 1.
set ytics 0, 2, yMax

set ls 10 lt 1 lc 0 lw 2.
set ls 11 lt 1 lc 9 lw 0.75

#~ tUnder='$\lambda=\frac{1}{10}\lambda_0^{}$'
#~ tOver='$\lambda=10\,\lambda_0^{}$'
#~ tGood_under='$\lambda_{\min}^{} + \text{twr}$'

tUnder = '{\small under-smeared\hspace{5.8em}}'
tOver = '{\small over-smeared\hspace{5.2em}}'
tGood_under='{\small natural resolution\hspace{7em}}'
tGood_over=tGood_under

set output 'Angular_3_smear_1deg.tex'
plot "./dat/A_ee_jjj_03_all_0.1_10.dat" using 1:4 ls 11 title tUnder , \
	'' using 1:3 ls 10 title tGood_under
	

set output 'Angular_3_smear_10deg.tex'
plot "./dat/A_ee_jjj_03_all_0.1_10.dat" using 1:3 ls 11 title tOver, \
	'' using 1:5 ls 10 title tGood_over

###############################

yMin = 1e-2
yMax = 2e2

set lmargin at screen 0.21

set ylabel offset 0.5

set format y '$10^{%L}$'
set mytics 10

set key at 0.85, 6e1

set logscale y
set yrange [1e-3:2e2]
set ytics 1e-2, 10, 1e3 add ('10' 10.) # , '1' 1.)

set output 'Angular_78_smear_1deg.tex'
plot "./dat/A_ee_jjj_78_all_0.1_10.dat" using 1:4 ls 11 title tUnder, \
	'' using 1:3 ls 10 title tGood_under

set output 'Angular_78_smear_10deg.tex'
plot "./dat/A_ee_jjj_78_all_0.1_10.dat" using 1:3 ls 11 title tOver, \
	'' using 1:5 ls 10 title tGood_over
	


