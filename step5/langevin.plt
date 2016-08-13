set term pngcairo
set out "langevin.png"

set xla "Time"
set style data line

set style data line

p "langevin.dat" t "Temperature" lw 3\
, "langevin.dat" u 1:3 t "Potential Energy" lw 3

