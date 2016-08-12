set term pngcairo
set out "energy.png"

set xla "Time"
set yla "Energy"
set key at 70,1.0
set style data line
p "energy.dat" u 1:2 lw 3 t "Kinetic Energy"\
, "energy.dat" u 1:3 lw 3 t "Potential Energy"\
, "energy.dat" u 1:4 lw 3 t "Total Energy"

