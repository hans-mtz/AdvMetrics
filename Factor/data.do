clear
set obs 10000
gen x0=1
gen x1=uniform()
gen z=rnormal()

gen em=rnormal()
gen e1=rnormal()*1.4
gen e0=rnormal()
gen ev=rnormal()

gen v1=1.0
gen v2=0.5
gen p1=0.5
gen p2=0.5
gen m1=-1.0
gen m2=-m1*p1/p2 // restrict the overall mean to be zero

gen theta=rnormal(m1,sqrt(v1))
gen u=uniform()
replace theta=rnormal(m2,sqrt(v2)) if u>p1
drop u

gen m=x0+x1+theta+em
gen y1=2*x0+x1+0.5*theta+e1
gen y0=x0-x1+2*theta+e0
gen v=x0-x1+z+theta+ev
gen d=0
replace d=1 if v>0
save data, replace
gen y=y1
replace y=y0 if d==0
outfile x0 x1 z m y d using data.raw, replace wide

