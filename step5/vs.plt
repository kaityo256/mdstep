set term pngcairo
set out "vs.png"

set xla "Time"
set style data line

set style data line

p "vs.dat" t "Temperature" lw 3\
, "vs.dat" u 1:3 t "Potential Energy" lw 3

