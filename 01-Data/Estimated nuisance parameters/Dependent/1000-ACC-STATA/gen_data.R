library(dplyr)
source("../helper_data_ACC.R")
set.seed(0); dat = data_mvn(1000)
dat$sim = 0
for (i in 1:1000) {
  set.seed(i)
  dat_temp = data_mvn(1000); dat_temp$sim = i
  dat = rbind(dat, dat_temp)
}
dat = dat %>%
  rename(R=D) %>%
  mutate(Xmiss = ifelse(R==1, X, NA)) %>%
  select(y,X,Xmiss,Z,R,sim)
write.csv(dat, "acc_example_stata.csv", row.names = FALSE, na = "")
