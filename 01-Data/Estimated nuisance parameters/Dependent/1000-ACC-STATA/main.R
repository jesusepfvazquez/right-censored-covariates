# Local simulation
rm(list = ls())

library(dplyr)
library(geex)
library(tictoc)
source("../helper_data_ACC.R")
source("../helper_mestimate_ACC.R")
dat = read.csv("acc_example_stata.csv")

simcalc = function(i=1){
  dat_temp = dat %>% subset(sim==1)
  
  est00 = estimate_beta(data_yXZ = dat_temp, model = y ~ Z+X)
  est10 = estimate_beta(data_yXZ = dat_temp, model = y ~ Z+W)
  est20 = estimate_beta_cc(data_yXZ = dat_temp, model = y ~ Z+W)
  
  # save results 
  myreturn = data.frame(cbind(c("oracle","naive","cc"),
                              rbind(est00$beta_est,est10$beta_est,est20$beta_est),
                              rbind(est00$se_est, est10$se_est, est20$se_est)))
  myreturn = data.frame(myreturn)
  colnames(myreturn) = c("Method",
                         "beta0.x", "beta1.x", "beta2.x","sigma.x", 
                         "beta0.y", "beta1.y", "beta2.y", "sigma.y")
  return(myreturn)
}


tic(); myresults = lapply(1:10, simcalc); toc()
# myresults = bind_rows(myresults)
myresults = bind_rows(myresults)
colnames(myresults) = c("Method","beta0.x", "beta1.x", "beta2.x", 
                             "beta0.y", "beta1.y", "beta2.y")

# Calculate bias 
theta = c(0,.2,.2)
myresults$beta0_coverage <- ((myresults$beta0.x + qnorm(0.975)*myresults$beta0.y) >= theta[1]) & 
  ((myresults$beta0.x - qnorm(0.975)*myresults$beta0.y) <= theta[1])
myresults$beta1_coverage <- ((myresults$beta1.x + qnorm(0.975)*myresults$beta1.y) >= theta[2]) & 
  ((myresults$beta1.x - qnorm(0.975)*myresults$beta1.y) <= theta[2])
myresults$beta2_coverage <- ((myresults$beta2.x + qnorm(0.975)*myresults$beta2.y) >= theta[3]) & 
  ((myresults$beta2.x - qnorm(0.975)*myresults$beta2.y) <= theta[3])

# est_threshold
est_threshold = 2

myresults %>%
  dplyr::select(Method,beta0.x, beta1.x, beta2.x, 
                beta0_coverage, beta1_coverage, beta2_coverage) %>%
  # calculate the means
  mutate(beta0_coverage = ifelse(beta0_coverage, 1, 0),
         beta1_coverage = ifelse(beta1_coverage, 1, 0),
         beta2_coverage = ifelse(beta2_coverage, 1, 0)) %>%
  group_by(Method) %>%
  summarise_all(mean) %>%
  print()

myresults %>%
  dplyr::select(Method,beta0.x, beta1.x, beta2.x) %>%
  group_by(Method) %>%
  summarise_all(var) %>%
  print()

