source("../../helper_data.R")
source("../../helper_mestimate_psi_V2.R")
args <- commandArgs(TRUE)
sim_seed <- 0 + as.integer(args)
set.seed(sim_seed)
n=1000

mybound = 10
var_c = 4
pcen=0.5
mu_c = qnorm(1-pcen)
print(n)

tic()

# generate data
mydata <- data_mvn(nSubjects = n,calculate_p = TRUE)
dat = mydata$mydata
myalpha = mydata$myalpha

est00 = estimate_beta(data_yXZ = dat, model = y ~ Z+AX)
est10 = estimate_beta(data_yXZ = dat, model = y ~ Z+AW)
est20 = estimate_beta_cc(data_yXZ = dat, model = y ~ Z+AW)
est30 = estimate_beta_ipw(data_yXZ = dat, model = y ~ Z+AW, myweight = "oracle")
est31 = estimate_beta_ipw(data_yXZ = dat, model = y ~ Z+AW, myweight = "mle")
est32 = estimate_beta_ipw_hat(data_yXZ = dat, model = y ~ Z+AW, myalpha)

est40 = estimate_beta_acc_lambda_alt(dat, y ~ Z+AW, myweight = "oracle")
est41 = estimate_beta_acc_lambda_hat(dat, y ~ Z+AW, myalpha = myalpha)

# est50=est00; est51=est00
est50 = estimate_beta_aipw_lambda_alt(dat, y ~ Z+AW, myweight = "oracle")
est51 = estimate_beta_aipw_lambda_hat(dat, y ~ Z+AW, myalpha = myalpha)

# save results
myreturn1 = data.frame(cbind(c("oracle","naive", "cc", "ipw", "ipw-hat", "ipw-hat2",
                               "acc-lambda", "acc-lambda-hat", 
                               "aipw-lambda", "aipw-lambda-hat"),
                             rbind(est00$beta_est,est10$beta_est,est20$beta_est,est30$beta_est, est31$beta_est, est32$beta_est,
                                   est40$beta_est, est41$beta_est,  
                                   est50$beta_est, est51$beta_est),
                             rbind(est00$se_est,est10$se_est,est20$se_est,est30$se_est, est31$se_est, est32$se_est,
                                   est40$se_est, est41$se_est,  
                                   est50$se_est, est51$se_est)
))
colnames(myreturn1) = c("Method",
                        "beta0.x", "beta1.x", "beta2.x","sigma.x",
                        "beta0.y", "beta1.y", "beta2.y")
                        
toc()

####################                        
# export beta estimates
write.csv(myreturn1, paste0("est_results_seed", sim_seed, ".csv"), row.names = F)

