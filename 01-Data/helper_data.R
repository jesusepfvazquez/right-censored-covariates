library(dplyr)
library(survival)
library(optimx)

# input 
## nSubjects: number of observations 
## dep_censoring: FALSE if cov(X,C|Z) = 0, and TRUE if != 0
## calculate_p: TRUE if we want pi(X,Z), else if we do not

data_mvn <- function(nSubjects = 10^3,  dep_censoring = TRUE){
  
  # defined parameters
  epsilon_sig = 1 #variance for the error component. (one value) 
  beta.vec = c(1,2,3) # true regression coefficients
  sig = 0.5 ## covariance between (X,Z) and (C,Z)
  
  ## generate subject IDs;
  b = 1:nSubjects
  
  ##############################################
  # Generate X, Z, C using multivariate normal #
  ##############################################
  mu_x = 0
  mu_z = 0
  mu_c=0
  var_x = 1; var_z = 1; var_c = 4
  var_xz = sig; var_cz = var_xz; var_xc =var_xz*var_cz
  if (dep_censoring) var_xc = 0.6
  
  # simulate using multivariate normal distribution
  Sigma= diag(3); diag(Sigma) = c(var_x, var_z, var_c);
  Mu = c(mu_x, mu_z, mu_c)
  Sigma[1,2] = var_xz; Sigma[1,3] = var_xc; Sigma[2,3] = var_cz
  Sigma[2,1] = var_xz; Sigma[3,1] = var_xc; Sigma[3,2] = var_cz
  
  # Generate data: X,Z,C
  dat = MASS::mvrnorm(n = nSubjects, mu = Mu, Sigma = Sigma,
                      tol = 1e-10, empirical = FALSE, EISPACK = FALSE)
  
  X = cbind(dat[,1])
  Z = cbind(dat[,2])
  C = cbind(dat[,3])
  A = rnorm(nSubjects, mean = 0, sd = 1)
  
  ##############################################################
  # Generate X, Z, C using conditional distribution based on Z #
  ##############################################################
  
  ## generate intercept
  int           = cbind(rep(1,nSubjects))
  
  ## simulate error terms
  er            = rnorm(nSubjects,mean=0,sd=epsilon_sig)
  
  # Calculate W, AW, AX, and delta, 
  W             = pmin(X, C)
  delta         = ifelse(X <= C,1,0)
  AX = A-X
  AW = A-W
  
  ## compute linear predictor space and outcomes
  x             = cbind(int,Z,AX) 
  linPred       = x%*%beta.vec
  
  ## y is a normal random variable 
  y             = linPred + er;
  
  # Create dataset 
  dat           = as.data.frame(cbind(y,Z,delta,X,C,W,b,A, AX, AW))
  colnames(dat) = c("y","Z","D","X","C","W","b","A", "AX", "AW")
  
  ##########################################################
  #################Calculate probabilities #################
  ##########################################################
  
  ###########################
  ## (X,Z)
  # oracle
  dat$myp_xz = pnorm(dat$X,
                     mean = mu_c + (var_x*var_z - var_xz^2)^(-1)*
                       ((var_xc*var_z-var_cz*var_xz)*(dat$X - mu_x) +
                          (var_cz*var_x-var_xc*var_xz)*(dat$Z - mu_z)),
                     sd = sqrt(var_c - (var_x*var_z - var_xz^2)^(-1)*
                                 ((var_xc*var_z-var_cz*var_xz)*(var_xc) +
                                    (var_cz*var_x-var_xc*var_xz)*(var_cz))),
                     lower.tail = FALSE)

  ##########################
  ## (Y,Z) 
  
  # oracle
  sdXZ = sqrt(var_x-var_xz^2/var_z)
  sdCZ = sqrt(var_c-var_cz^2/var_z)
  
  myp_yz.b = function(data.){
    
    # data. = dat[1,]
    meanXZ = mu_x + (var_xz/var_z)*(data.$Z-mu_z)
    integral_func_num.b = function(t){
      val = rep(NA, length(t))
      for (i in 1:length(t)){
        myp = pnorm(t[i],
                    mean = mu_c + (var_x*var_z - var_xz^2)^(-1)*
                      ((var_xc*var_z-var_cz*var_xz)*(t[i] - mu_x) +
                         (var_cz*var_x-var_xc*var_xz)*(data.$Z - mu_z)),
                    sd = sqrt(var_c - (var_x*var_z - var_xz^2)^(-1)*
                                ((var_xc*var_z-var_cz*var_xz)*(var_xc) +
                                   (var_cz*var_x-var_xc*var_xz)*(var_cz))),
                    lower.tail = FALSE)
        val[i] = myp*
          dnorm(data.$y -  cbind(1,data.$Z,data.$A - t[i])%*%beta.vec, 0, sd = epsilon_sig)*
          dnorm(t[i], mean = meanXZ, sd = sdXZ)
      }
      return(val)
    }
    myp_yz_num = integrate(integral_func_num.b, lower = -Inf, upper = Inf, 
                       rel.tol = 1e-7, abs.tol = 1e-7)$value
    
    integral_func_denom.b = function(t){
      val = rep(NA, length(t))
      for (i in 1:length(t)){
        val[i] = dnorm(data.$y -  cbind(1,data.$Z,data.$A - t[i])%*%beta.vec, 0, sd = epsilon_sig)*
          dnorm(t[i], mean = meanXZ, sd = sdXZ)
      }
      return(val)
    }
    myp_yz_denom = integrate(integral_func_denom.b, lower = -Inf, upper = Inf, 
                           rel.tol = 1e-7, abs.tol = 1e-7)$value
    return(myp_yz_num/myp_yz_denom)
  } 
  dat$myp_yz = lapply(1:nrow(dat), function(t.) myp_yz.b(dat[t.,])) %>% unlist()

  
  return(dat)
}

# data generating function used to generate data from Bartlette 2024.
data_mvn_bartlette <- function(nSubjects = 10^3){
  
  ## parameter values from function
  # nSubjects = 1000
  pcen = 0.5 # percent right-censored
    
  ## generate subject IDs;
  b = 1:nSubjects
  
  # Step 1: simulate delta
  delta = rbinom(nSubjects, 1, pcen)
  
  # Step 2: Simulate (X,Z,Y)|delta
  gamma_X0 = 0; gamma_Z0 = 0; gamma_Y0 = 0
  gamma_X1 = 1; gamma_Z1 = 1;
  var_x = 1; var_z = 1; var_y = 1
  var_xz = 0.25; var_zy = 0.25; var_xy = 0.25
  
  # condition for delta \perp Y | X, Z
  gamma_Y1 = (gamma_X1*(var_xy*var_z^2 - var_xz*var_zy) + 
                gamma_Z1*(var_zy*var_x^2 - var_xy*var_xz))/(var_x^2*var_z^2 - var_xz^2)
  
  # assigning mu
  mu_x = gamma_X0 +gamma_X1*delta
  mu_z = gamma_Z0 +gamma_Z1*delta
  mu_y = gamma_Y0 +gamma_Y1*delta
  
  # simulate using multivariate normal distribution
  Mu = cbind(0,0,0)
  Sigma= diag(3); diag(Sigma) = c(var_x, var_z, var_y);
  Sigma[1,2] = var_xz; Sigma[1,3] = var_xy; Sigma[2,3] = var_zy
  Sigma[2,1] = var_xz; Sigma[3,1] = var_xy; Sigma[3,2] = var_zy
  
  # Generate data: X,Z,Y
  dat = MASS::mvrnorm(n = nSubjects, mu = Mu, Sigma = Sigma,
                      tol = 1e-10, empirical = FALSE, EISPACK = FALSE)  + cbind(mu_x,mu_z,mu_y)
  
  # Step 4: generate C such that C \!perp! X | Z
  C = dat[,1] - rbeta(nSubjects, exp(dat[,2]), 1)
  W = ifelse(delta == 1, dat[,1], C)
  
  # save results and return 
  dat = data.frame(b, dat, delta, W)
  colnames(dat) = c("b", "X", "Z", "y", "D", "W")
  
  # calculate oracle weights 
  dat$myp = glm(D ~ X+Z, data = dat, family = binomial(link= "logit"))$fitted
  dat$myp_yz = glm(D ~ y+Z, data = dat, family = binomial(link= "logit"))$fitted
  
  # calculate random weights
  dat$myp_unif = runif(nSubjects, min = 0.1, max = 1)
  
  ## beta values: Y | X + Z
  beta1 = (var_xy*var_z^2 - var_zy*var_xz)/(var_z^2*var_x^2-var_xz^2)
  beta2 = (var_zy*var_x^2 - var_xy*var_xz)/(var_x^2*var_z^2-var_xz^2)
  beta0 = gamma_Y1 - beta1*gamma_X1 - beta2*gamma_Z1
  theta1 = c(beta0, beta1, beta2)
  # print(theta1)
  
  ## gamma values: X|Z,Y,delta=1
  dz = (var_xy*var_z^2 - var_xz*var_zy)/(var_z^2*var_y^2 - var_zy^2)
  dy = (var_xz*var_y^2 - var_xy*var_zy)/(var_y^2*var_z^2 - var_zy^2)
  d0 = gamma_X0 + gamma_X1 - dz*(gamma_Z0 + gamma_Z1) - dy*(gamma_Y0 - gamma_Y1)
  psi2 = var_x^2 - (dz*var_xz + dy*var_xy)
  theta2 = c(d0,dz,dy,psi2)
  # print(theta2)
  
  dat$meanX = theta2[1] + theta2[2]*dat$Z + theta2[3]*dat$y
  dat$meanX_wrong  = 0
  
  return(dat)
}


