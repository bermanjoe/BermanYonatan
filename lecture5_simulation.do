clear

* Setting the number of simulated values
set obs 10000

* Running variable - randomly drawn from a normal distribution aroudn 50
* X is some test score
gen X = rnormal(50, 25)

* We don't want negative test scores, or higher than 100, we assume these are errors and drop them
drop if X < 0
drop if X > 100

* Generate treatment vector
gen T = 0
replace T = 1 if X > 50

* Generate errors
gen errs = rnormal(0, 20)

* Generate values
gen y = 25 + 1.5*X + 40 * T + errs

* Plot values
scatter y X if T==0, msize(vsmall) || scatter y X if T==1, msize(vsmall) legend(off) xline(50, lstyle(foreground)) || lfit y X if T ==0, color(red) || lfit y X if T ==1, color(red) ytitle("Outcome (Y)")  xtitle("Test Score (X)") 

* Create transformed x (i.e., x-c0)
gen tx = X - 50

gen txT = tx * T

* Estimate treatment effect
reg y tx txT T

* regress without transformation
gen XT = X * T
reg y tx txT T