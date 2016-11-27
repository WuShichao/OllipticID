

f(x)=a+b*x

set fit logfile "CAfit.log"

set fit errorvariables

fit f(x) "vL1.dat" via a,b



set ylabel "Log(E(h)) (adim)"

set xlabel "Log(h) (adim)"

set key bottom


set title sprintf("E(h)~C h^ p :\n Log(C)= %f +/- %f \n p = %f +/- %f",a,a_err,b,b_err)

plot "vL1.dat" t "Olliptic", a+b*x t sprintf("%f + %f x",a,b)

pause -1

