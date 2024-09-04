source("../../helper_data.R")
source("../../helper_mestimate_psi_V3.R")
args <- commandArgs(TRUE)
sim_seed <- 0 + as.integer(args)
set.seed(sim_seed)
n=300
pcen = 0.5
mysig = 0.5
mybound = 10
var_c = 4
mu_c = qnorm(1-pcen)
dep_censoring. = FALSE
print(n)

tic()

# generate data
dat = data_mvn_c(nSubjects = n, sig = mysig,  dep_censoring = dep_censoring.)
print(mean(dat$D))

est00 = estimate_beta(data_yXZ = dat, model = y ~ Z+AX)
est10 = estimate_beta(data_yXZ = dat, model = y ~ Z+AW)
est20 = estimate_beta_cc(data_yXZ = dat, model = y ~ Z+AW)

est30 = estimate_beta_ipw(data_yXZ = dat, model = y ~ Z+AW, myweight = "oracle")
est31 = estimate_beta_ipw(data_yXZ = dat, model = y ~ Z+AW, myweight = "uniform")
est32 = estimate_beta_ipw(data_yXZ = dat, model = y ~ Z+AW, myweight = "yz")

# est20 = est32
# est00 = est20; est10 = est20; est30 = est20; est31 = est20
# est40 = est00; est41 = est00; est42 = est00; est43=est00; est44 = est00; est45=est00;est46=est00;est47=est00
# est50 = est00; est51 = est00; est52 = est00; est53=est00; est54 = est00; est55=est00; est56=est00; est57=est00
# est60 = est00; est61 = est00; est62 = est00; est63=est00; est64 = est00; est65=est00; est66=est00; est67=est00
#est70 = est00; est71 = est00

print("starting acc")

# acc
est40 = estimate_beta_acc(dat, y ~ Z+AW, eta1= FALSE, myweight = "oracle", dep_censoring = dep_censoring.)
est41 = estimate_beta_acc(dat, y ~ Z+AW, eta1= TRUE, myweight = "oracle", dep_censoring = dep_censoring.)
est42 = estimate_beta_acc(dat, y ~ Z+AW, eta1= FALSE, myweight = "uniform", dep_censoring = dep_censoring.)
est43 = estimate_beta_acc(dat, y ~ Z+AW, eta1= TRUE, myweight = "uniform", dep_censoring = dep_censoring.)
est44 = estimate_beta_acc_lambda_close(dat, y ~ Z+AW, myweight = "oracle")
est45 = estimate_beta_acc_lambda_close(dat, y ~ Z+AW, myweight = "uniform")

print("starting macc")

# macc
est50 = estimate_beta_macc(dat, y ~ Z+AW, eta1= FALSE, myweight = "oracle", dep_censoring = dep_censoring.)
est51 = estimate_beta_macc(dat, y ~ Z+AW, eta1= TRUE, myweight = "oracle", dep_censoring = dep_censoring.)
est52 = estimate_beta_macc(dat, y ~ Z+AW, eta1= FALSE, myweight = "uniform", dep_censoring = dep_censoring.)
est53 = estimate_beta_macc(dat, y ~ Z+AW, eta1= TRUE, myweight = "uniform", dep_censoring = dep_censoring.)
est54 = estimate_beta_macc_lambda_close(dat, y ~ Z+AW, myweight = "oracle")
est55 = estimate_beta_macc_lambda_close(dat, y ~ Z+AW, myweight = "uniform")

print("starting aipw")

# aipw
est60 = estimate_beta_aipw(dat, y ~ Z+AW, eta1= FALSE, myweight = "oracle", dep_censoring = dep_censoring.)
est61 = estimate_beta_aipw(dat, y ~ Z+AW, eta1= TRUE, myweight = "oracle", dep_censoring = dep_censoring.)
est62 = estimate_beta_aipw(dat, y ~ Z+AW, eta1= FALSE, myweight = "uniform", dep_censoring = dep_censoring.)
est63 = estimate_beta_aipw(dat, y ~ Z+AW, eta1= TRUE, myweight = "uniform", dep_censoring = dep_censoring.)
est64 = estimate_beta_aipw_lambda_close(dat, y ~ Z+AW, myweight = "oracle")
est65 = estimate_beta_aipw_lambda_close(dat, y ~ Z+AW, myweight = "uniform")


# save results
myreturn1 = data.frame(cbind(c("oracle","naive", "cc", "ipw", "ipw-w", "ipw-yz", #"ipw-yz-w",
                              "acc", "acc-eta1", "acc-w", "acc-eta1-w", "acc-lambda", "acc-lambda-w", 
                              "macc", "macc-eta1", "macc-w", "macc-eta1-w", "macc-lambda", "macc-lambda-w", 
                              "aipw", "aipw-eta1", "aipw-w", "aipw-eta1-w", "aipw-lambda", "aipw-lambda-w"),
                            rbind(est00$beta_est,est10$beta_est,est20$beta_est,est30$beta_est, est31$beta_est, est32$beta_est, #est33$beta_est,
                                  est40$beta_est, est41$beta_est, est42$beta_est, est43$beta_est, est44$beta_est, est45$beta_est, 
                                  est50$beta_est, est51$beta_est, est52$beta_est, est53$beta_est, est54$beta_est, est55$beta_est,
                                  est60$beta_est, est61$beta_est, est62$beta_est, est63$beta_est, est64$beta_est, est65$beta_est),
                            rbind(est00$se_est,est10$se_est,est20$se_est,est30$se_est,est31$se_est,est32$se_est, #est33$se_est,
                                  est40$se_est, est41$se_est, est42$se_est, est43$se_est, est44$se_est, est45$se_est, 
                                  est50$se_est, est51$se_est, est52$se_est, est53$se_est, est54$se_est, est55$se_est,
                                  est60$se_est, est61$se_est, est62$se_est, est63$se_est, est64$se_est, est65$se_est)
))
colnames(myreturn1) = c("Method",
                        "beta0.x", "beta1.x", "beta2.x","sigma.x",
                        "beta0.y", "beta1.y", "beta2.y")

print("starting mle")

#mle

est70 = estimate_beta_likelihood(dat, y ~ Z+AW, eta1= FALSE, dep_censoring = dep_censoring.)
#est71 = est70
est71 = estimate_beta_likelihood(dat, y ~ Z+AW, eta1= TRUE, dep_censoring = dep_censoring.)

myreturn2 = data.frame(cbind(c("mle", "mle-w"),
                            rbind(est70$beta_est,est71$beta_est),
                            rbind(est70$se_est,est71$se_est)))[,-9]
colnames(myreturn2) = c("Method",
                       "beta0.x", "beta1.x", "beta2.x","sigma.x",
                       "beta0.y", "beta1.y", "beta2.y")

# save results
myreturn = rbind(myreturn1, myreturn2)
#myreturn = myreturn1
toc()

####################                        
# export beta estimates
#myreturn = rbind(myreturn1, myreturn2)
write.csv(myreturn, paste0("est_results_seed", sim_seed, ".csv"), row.names = F)

