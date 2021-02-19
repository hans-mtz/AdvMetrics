clear
set mem 100m
set obs 50253
gen x1=1
gen x2=uniform()*10
gen e=invnorm(uniform())
egen me=mean(e)
replace e=e-me
replace e=e*2
gen y=-x1+2.36*x2+e

reg y x1 x2

outfile y x1 x2 e using data, replace w
