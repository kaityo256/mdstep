set term png 18
set out "temperature.png"
set xla "Time"
set yla "Temperature"
unset key
set xtics 10
set ytics 0.2
set style data line
p [][0:]"output.dat" u 1:2 lw 3
set out "pressure.png"
set yla "Pressure"
set ytics 0.2
p [][-1:0.2]"output.dat" u 1:3 lw 3

