set term pngcairo
set out "nosehoover.png"

set xla "Time"
set style data line

set style data line

set ytics 0.5
p "nosehoover.dat" t "Temperature" lw 3\
, "nosehoover.dat" u 1:3 t "Potential Energy" lw 3

