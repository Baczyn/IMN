set terminal post enhanced colour solid font 20



unset key
set xlabel "x"
set ylabel "y"
set grid

set output "RK2.png"
plot "RK2.txt" u 1:2 w lp,"RK2.txt" u 1:3 w lp

set output "picard.png"
plot "picard.txt" u 1:2 w lp,"picard.txt" u 1:3 w lp

set output "newton.png"
plot "newton.txt" u 1:2 w lp,"newton.txt" u 1:3 w lp

