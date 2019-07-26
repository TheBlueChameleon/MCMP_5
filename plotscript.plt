# =========================================================================== #
# global output properties

set terminal pdf
set output "results.pdf"

# =========================================================================== #
# calibration plots

set dgrid3d 20,50 
set xyplane at 0

set logscale x 2

set palette rgbformulae 23,28,3
unset surface
set pm3d at b
set pm3d 
set pm3d map

set xrange [2:16]

set xlabel 'leapfrog steps n'
set ylabel 'leapfrog time step dt'

show title
set key off

# --------------------------------------------------------------------------- #
set title "Computation Cost for low T\nEstimate for t_{indep}\n"
splot './out/cal.dat' u 1:2:3

set title "Computation Cost for high T\nEstimate for t_{indep}"
splot './out/cal.dat' u 4:5:6

# =========================================================================== #
# behaviour plots

reset
set nokey
show title

set  xlabel 'Temperature T'
show xlabel

set style line 1 linecolor rgb "red"   pointtype 1   # plus
set style line 2 linecolor rgb 'blue'  pointtype 1 
set style line 3 linecolor rgb 'green' pointtype 1

# --------------------------------------------------------------------------- #
set title 'Energy Density vs. Temperature'

set  ylabel 'e(T)'
show ylabel

plot  './out/EMCX.dat'      using 1:2:3 title 'data'      with errorbars linestyle 1

# --------------------------------------------------------------------------- #
set title 'Magnetization Density vs. Temperature'

set  ylabel '|m(T)|'
show ylabel

plot  './out/EMCX.dat'      using 1:4:5 title 'data'      with errorbars linestyle 1

# --------------------------------------------------------------------------- #

set title 'Specific Heat Density vs. Temperature'

set  ylabel 'c(T)'
show ylabel

plot  './out/EMCX.dat'      using 1:6:7 title 'data'      with errorbars linestyle 1

# --------------------------------------------------------------------------- #

set title 'Magnetic Susceptibility Density vs. Temperature'

set  ylabel '{/Symbol c}(T)'
show ylabel

plot  './out/EMCX.dat'      using 1:8:9 title 'data'      with errorbars linestyle 1

# --------------------------------------------------------------------------- #
