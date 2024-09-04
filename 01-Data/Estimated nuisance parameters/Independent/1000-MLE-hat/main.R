source("../../helper_data.R")
source("../../helper_mestimate_psi_V2.R")
args <- commandArgs(TRUE)
sim_seed <- 0 + as.integer(args)
set.seed(sim_seed)
n=1000
pcen =0.5
mysig = 0.5
mybound = 10
var_c = 4
mu_c = qnorm(1-pcen)
print(n)

tic()

# generate data
mydata = data_mvn(nSubjects = n, sig = mysig,  dep_censoring = FALSE, calculate_p = TRUE)
dat = mydata$mydata
myalpha = mydata$myalpha

est00 = estimate_beta(data_yXZ = dat, model = y ~ Z+AX)
est10 = estimate_beta(data_yXZ = dat, model = y ~ Z+AW)
est20 = estimate_beta_cc(data_yXZ = dat, model = y ~ Z+AW)
est30 = estimate_beta_ipw(data_yXZ = dat, model = y ~ Z+AW, myweight = "oracle")
est31 = estimate_beta_ipw(data_yXZ = dat, model = y ~ Z+AW, myweight = "uniform")

#print("starting acc")

#est40 = est00; est41 = est00; est42 = est00; est43=est00; est44 = est00; est45=est00; est46 = est00; est47=est00
#est50 = est00; est51 = est00; est52 = est00; est53=est00; est54 = est00; est55=est00; est56 = est00; est57=est00
#est60 = est00; est61 = est00

# acc
#est40 = estimate_beta_acc(dat, y ~ Z+AW, eta1= FALSE, myweight = "oracle", dep_censoring = FALSE)
#est41 = estimate_beta_acc(dat, y ~ Z+AW, eta1= TRUE, myweight = "oracle", dep_censoring = FALSE)
#est42 = estimate_beta_acc(dat, y ~ Z+AW, eta1= FALSE, myweight = "uniform", dep_censoring = FALSE)
#est43 = estimate_beta_acc(dat, y ~ Z+AW, eta1= TRUE, myweight = "uniform", dep_censoring = FALSE)
#est44 = estimate_beta_acc_lambda(dat, y ~ Z+AW, eta1= FALSE, myweight = "oracle", dep_censoring = FALSE)
#est45 = estimate_beta_acc_lambda_alt(dat, y ~ Z+AW, myweight = "oracle")
#est46 = estimate_beta_acc_lambda(dat, y ~ Z+AW, eta1= FALSE, myweight = "uniform", dep_censoring = FALSE)
#est47 = estimate_beta_acc_lambda_alt(dat, y ~ Z+AW, myweight = "uniform")

#print("starting aipw")

#aipw
#est50 = estimate_beta_aipw(dat, y ~ Z+AW, eta1= FALSE, myweight = "oracle", dep_censoring = FALSE)
#est51 = estimate_beta_aipw(dat, y ~ Z+AW, eta1= TRUE, myweight = "oracle", dep_censoring = FALSE)
#est52 = estimate_beta_aipw(dat, y ~ Z+AW, eta1= FALSE, myweight = "uniform", dep_censoring = FALSE)
#est53 = estimate_beta_aipw(dat, y ~ Z+AW, eta1= TRUE, myweight = "uniform", dep_censoring = FALSE)
#est54 = estimate_beta_aipw_lambda(dat, y ~ Z+AW, eta1= FALSE, myweight = "oracle", dep_censoring = FALSE)
#est55 = estimate_beta_aipw_lambda_alt(dat, y ~ Z+AW, myweight = "oracle")
#est56 = estimate_beta_aipw_lambda(dat, y ~ Z+AW, eta1= FALSE, myweight = "uniform", dep_censoring = FALSE)
#est57 = estimate_beta_aipw_lambda_alt(dat, y ~ Z+AW, myweight = "uniform")

#myreturn1 = data.frame(cbind(c("oracle"),
#				rbind(est00$beta_est

# save results
  myreturn1 = data.frame(cbind(c("oracle","naive", "cc", "ipw", "ipw-w"),
#                                "acc", "acc-eta1", "acc-w", "acc-eta1-w", "acc-lambda", "acc-lambda-eta1", "acc-lambda-w", "acc-lambda-eta1-w",
#                                "aipw", "aipw-eta1", "aipw-w", "aipw-eta1-w", "aipw-lambda", "aipw-lambda-eta1",  "aipw-lambda-w", "aipw-lambda-eta1-w"),
                              rbind(est00$beta_est,est10$beta_est,est20$beta_est,est30$beta_est, est31$beta_est),
#                                    est40$beta_est, est41$beta_est, est42$beta_est, est43$beta_est, est44$beta_est, est45$beta_est, est46$beta_est, est47$beta_est,
#                                    est50$beta_est, est51$beta_est, est52$beta_est, est53$beta_est, est54$beta_est, est55$beta_est, est56$beta_est, est57$beta_est),
                              rbind(est00$se_est,est10$se_est,est20$se_est,est30$se_est,est31$se_est)
#                                    est40$se_est, est41$se_est, est42$se_est, est43$se_est, est44$se_est, est45$se_est, est46$se_est, est47$se_est,
#                                    est50$se_est, est51$se_est, est52$se_est, est53$se_est, est54$se_est, est55$se_est, est56$se_est, est57$se_est)
  ))
  colnames(myreturn1) = c("Method",
                          "beta0.x", "beta1.x", "beta2.x","sigma.x",
                          "beta0.y", "beta1.y", "beta2.y")

####################  
print("starting mle")

# mle
est60 = estimate_beta_likelihood_hat(dat, y ~ Z+AW, myalpha = myalpha)
est61=est60
#est61 = estimate_beta_likelihood(dat, y ~ Z+AW, eta1= TRUE, dep_censoring = FALSE)

# save results
 myreturn2 = data.frame(cbind(c("mle", "mle-w"),
                              rbind(est60$beta_est,est61$beta_est),
                              rbind(est60$se_est,est61$se_est)))[,-9]
  colnames(myreturn2) = c("Method",
                         "beta0.x", "beta1.x", "beta2.x","sigma.x",
                         "beta0.y", "beta1.y", "beta2.y")

toc()

####################                        
# export beta estimates
myreturn = rbind(myreturn1, myreturn2)
write.csv(myreturn, paste0("est_results_seed", sim_seed, ".csv"), row.names = F)

print(est20)
