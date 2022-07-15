clear all
qui set obs 100000
gen x = 2*3.14159265*runiform()
gen y = sin(x) + rnormal()/10
forval i = 1/10 {
    timer on `i'
    lpoly y x, n(1000) kernel(epan2) degree(1) nograph
    timer off `i'
}
timer list
