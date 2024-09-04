// jesus vazquez
// 04/18/2024

clear
import delimited "acc_example_stata.csv"

keep if sim == `1'
gen mysim`1' = 0

// gen imp1 = rnormal(0,1)
// gen imp2 = rnormal(0,1)
// gen imp3 = rnormal(0,1)
// gen imp4 = rnormal(0,1)
// gen imp5 = rnormal(0,1)
// gen imp6 = rnormal(0,1)
// gen imp7 = rnormal(0,1)
// gen imp8 = rnormal(0,1)
// gen imp9 = rnormal(0,1)
// gen imp10 = rnormal(0,1)
// gen imp11 = rnormal(0,1)
// gen imp12 = rnormal(0,1)
// gen imp13 = rnormal(0,1)
// gen imp14 = rnormal(0,1)
// gen imp15 = rnormal(0,1)
// gen imp16 = rnormal(0,1)
// gen imp17 = rnormal(0,1)
// gen imp18 = rnormal(0,1)
// gen imp19 = rnormal(0,1)
// gen imp20 = rnormal(0,1)
// gen imp21 = rnormal(0,1)
// gen imp22 = rnormal(0,1)
// gen imp23 = rnormal(0,1)
// gen imp24 = rnormal(0,1)
// gen imp25 = rnormal(0,1)
// gen imp26 = rnormal(0,1)
// gen imp27 = rnormal(0,1)
// gen imp28 = rnormal(0,1)
// gen imp29 = rnormal(0,1)
// gen imp30 = rnormal(0,1)
//
// // correct weight 
// augcca z, y(y) x(xmiss) imps(imp1-imp30) pivars(y z) mysim(mysim`1')
//
// // incorrect weight specification
// augcca z, y(y) x(xmiss) imps(imp1-imp30) pivars(z) mysim(mysim`1')

// incorrect augmentation
gen imp31 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp32 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp33 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp34 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp35 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp36 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp37 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp38 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp39 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp40 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp41 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp42 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp43 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp44 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp45 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp46 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp47 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp48 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp49 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp50 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp51 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp52 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp53 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp54 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp55 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp56 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp57 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp58 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp59 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y
gen imp60 = rnormal(0,0.9) + 0.88 + 0.20*z + 0.20*y

// correct weight 
augcca z, y(y) x(xmiss) imps(imp31-imp60) pivars(y z) mysim(mysim`1')

// incorrect weight specification
augcca z, y(y) x(xmiss) imps(imp31-imp60) pivars(z) mysim(mysim`1')




