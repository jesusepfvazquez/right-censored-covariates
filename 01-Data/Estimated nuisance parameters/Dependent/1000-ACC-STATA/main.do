// jesus vazquez
// 04/18/2024

clear
import delimited "acc_example_stata.csv"

keep if sim == `1'
gen mysim`1' = 0

gen imp1 = rnormal(0,1)
gen imp2 = rnormal(0,1)
gen imp3 = rnormal(0,1)
gen imp4 = rnormal(0,1)
gen imp5 = rnormal(0,1)
gen imp6 = rnormal(0,1)
gen imp7 = rnormal(0,1)
gen imp8 = rnormal(0,1)
gen imp9 = rnormal(0,1)
gen imp10 = rnormal(0,1)
gen imp11 = rnormal(0,1)
gen imp12 = rnormal(0,1)
gen imp13 = rnormal(0,1)
gen imp14 = rnormal(0,1)
gen imp15 = rnormal(0,1)
gen imp16 = rnormal(0,1)
gen imp17 = rnormal(0,1)
gen imp18 = rnormal(0,1)
gen imp19 = rnormal(0,1)
gen imp20 = rnormal(0,1)
gen imp21 = rnormal(0,1)
gen imp22 = rnormal(0,1)
gen imp23 = rnormal(0,1)
gen imp24 = rnormal(0,1)
gen imp25 = rnormal(0,1)
gen imp26 = rnormal(0,1)
gen imp27 = rnormal(0,1)
gen imp28 = rnormal(0,1)
gen imp29 = rnormal(0,1)
gen imp30 = rnormal(0,1)

// correct weight 
augcca z, y(y) x(xmiss) imps(imp1-imp30) pivars(y z) mysim(mysim`1')

// incorrect weight specification
augcca z, y(y) x(xmiss) imps(imp1-imp30) pivars(z) mysim(mysim`1')


