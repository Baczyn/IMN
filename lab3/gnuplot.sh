set terminal post enhanced colour solid font 20


unset key
set xlabel "x"
set ylabel "y"
set grid

set output "RK2_1.eps"
plot "RK2.txt" u 1:2 w lp

set output "RK2_2.eps"
plot "RK2.txt" u 1:3 w lp

set output "RK2_3.eps"
plot "RK2.txt" u 1:3 w lp

set output "RK2_4.eps"
plot "RK2.txt" u 3:4 w lp

set output "trapez_1.eps"
plot "trapez.txt" u 1:2 w lp

set output "trapez_2.eps"
plot "trapez.txt" u 1:3 w lp

set output "trapez_3.eps"
plot "trapez.txt" u 1:4 w lp

set output "trapez_4.eps"
plot "trapez.txt" u 3:4 w lp
