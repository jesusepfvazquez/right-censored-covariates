## m-estimating equations & other helper functions ##

#############################
# known nuisance parameters #
#############################

# oracle or naive
estimate_beta <- function(data_yXZ, model){
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    # data = dat[1,]
    # X <- model.matrix(object = y ~ AW + Z, data = data)
    # Y <- model.response(model.frame(formula = y ~ AW + Z, data = data))
    # theta = 1:3
    # psi=1
    
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    
    function(theta){
      
      # first component of score
      myscore <- function(){
        # useful quantities
        mu <- X %*% theta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        return(c(dotbeta))
      }
      
      return(myscore())
    }
  }
  
  ######## EST INITIAL PARAMETER ESTIMATES ##########
  mylmer  = lm(model, data = data_yXZ, weight = D)
  g = c(coef(mylmer))
  psi = summary(mylmer)$sigma
  # print(g); print(psi)
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    root_control = setup_root_control(start = c(g)),
    outer_args = list(formula = model),
    deriv_control = setup_deriv_control(method = 'Richardson')
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% matrix()
  # print(beta_estimates)
  # head(model.matrix(object = model, data = data_yXZ))

  psi_updated = sum((data_yXZ$y - model.matrix(object = model, data = data_yXZ)%*%beta_estimates)^2/(nrow(data_yXZ))) %>% sqrt
  beta_estimates = rbind(beta_estimates, psi_updated)
  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% matrix()
  return(list(beta_est = t(beta_estimates), se_est = t(se_estimates)))
}

# complete case
estimate_beta_cc <-  function(data_yXZ, model){
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    D <- data$D
    
    function(theta){
      
      # first component of score
      myscore <- function(){
        # useful quantities
        mu <- X %*% theta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        return(D*c(dotbeta))
      }
      
      return(myscore())
    }
  }
  
  ######## EST INITIAL PARAMETER ESTIMATES ##########
  mylmer  = lm(model, data = data_yXZ, weight = D)
  g = c(coef(mylmer))
  psi = summary(mylmer)$sigma
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    root_control = setup_root_control(start = c(g)),
    outer_args = list(formula = model),
    deriv_control = setup_deriv_control(method = 'Richardson')
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% matrix()
  data_yXZ = data_yXZ %>% subset(D==1)  
  psi_updated = sum((data_yXZ$y - model.matrix(object = model, data = data_yXZ)%*%beta_estimates)^2/(nrow(data_yXZ))) %>% sqrt
  beta_estimates = rbind(beta_estimates, psi_updated)
  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% matrix()
  return(list(beta_est = t(beta_estimates), se_est = t(se_estimates)))
}

# ipw
estimate_beta_ipw <- function(data_yXZ, model, myweight = "oracle"){
  
  #######################################################
  # define weight components #
  #######################################################
  data_yXZ$myp = data_yXZ$myp_xz
  if (myweight  == "cox_xz") data_yXZ$myp = data_yXZ$myp_xz_cox
  if (myweight  == "mle") data_yXZ$myp = data_yXZ$myp_xz_mle
  if (myweight  == "uniform") data_yXZ$myp = data_yXZ$myp_uniform
  if (myweight  == "yz") data_yXZ$myp = data_yXZ$myp_yz
  
  #######################################################
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    D <- data$D
    myp <- data$myp
    
    function(theta){
      
      # first component of score
      CCscore <- function(){
        # useful quantities
        mu <- X %*% theta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        return(c(dotbeta))
      }
      
      return(CCscore()*D/myp)
    }
  }
  
  ######## EST INITIAL PARAMETER ESTIMATES ##########
  mylmer  = lm(model, data = data_yXZ, weights = D/myp)
  g = c(coef(mylmer))
  psi = summary(mylmer)$sigma
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    root_control = setup_root_control(start = c(g)),
    outer_args = list(formula = model),
    deriv_control = setup_deriv_control(method = 'Richardson')
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% matrix()
  data_yXZ = data_yXZ %>% subset(D==1)  
  psi_updated = sum((data_yXZ$y - model.matrix(object = model, data = data_yXZ)%*%beta_estimates)^2/(nrow(data_yXZ))) %>% sqrt
  beta_estimates = rbind(beta_estimates, psi_updated)
  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% matrix()
  return(list(beta_est = t(beta_estimates), se_est = t(se_estimates)))
}

# acc
estimate_beta_acc <- function(data_yXZ, model, eta1 = FALSE, myweight = "oracle", dep_censoring = FALSE){
  
  # data_yXZ = dat
  ######################################################
  # define parameters values for integration component #
  ######################################################
  sig = 0.5
  mu_x = 0; mu_z = 0; mu_c =0
  var_x = 1; var_z = 1; var_c = 1
  var_xz = sig; var_cz = var_xz; var_xc = var_xz*var_cz 
  sdXZ = sqrt(var_x-var_xz^2/var_z)
  sdCXZ = sqrt(var_c-var_cz^2/var_z)
  
  if (dep_censoring){
    var_xc = 0.6
    sdCXZ = sqrt(var_c - (var_x*var_z - var_xz^2)^(-1)*
                   ((var_xc*var_z-var_cz*var_xz)*(var_xc) +
                      (var_cz*var_x-var_xc*var_xz)*(var_cz)))
  }
  
  ################################
  # define eta1 misspecification #
  ################################
  mymeanXZ = function(data.){mu_x + (var_xz/var_z)*(data.$Z-mu_z)}
  
  if (myweight  == "oracle") data_yXZ$myp = data_yXZ$myp_yz
  if (myweight  == "uniform") data_yXZ$myp = data_yXZ$myp_uniform
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    W <- model.matrix(object = ~ A + Z, data = data)
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    D <- data$D
    z <- data$Z
    myp <- data$myp
    
    # need these for integration components
    meanXZ = mymeanXZ(data)
    
    function(theta){
      
      #################
      # 1st Component #
      #################
      CCscore <- function(){
        # useful quantities
        mu <- X %*% theta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        return(c(dotbeta))
      }
      
      #######################
      # Augmented Component #
      #######################
      # augmentation component
      psi_hat_i <- function(){
        
        # top integral
        likelihood_int = function(t=1){
          return(dnorm(Y -  cbind(W[,1],W[,2]-t, W[,3]) %*% theta, 0, sd = psi))
        }
        
        # score of integral
        score_int <- function(t=1,j=1){
          # generate empty array of values to save
          # create temporary X matrix
          Wtemp <- cbind(W[,1],W[,2]-t, W[,3])
          mu <- Wtemp %*% theta
          e = Y - mu
          # score equations for theta
          dotbeta = (e/psi^2)%*%Wtemp
          return(c(dotbeta)[j])
        }
        
        #top integral#
        integral_func_num <- function(t=1,j=1){
          val = rep(NA, length(t))
          meanCXZ = mu_c + (var_x*var_z - var_xz^2)^(-1)*
            ((var_xc*var_z-var_cz*var_xz)*(t - mu_x) +
               (var_cz*var_x-var_xc*var_xz)*(z - mu_z))
          for (i in 1:length(t)){
            val[i] = score_int(t[i],j=j)*
              likelihood_int(t[i])* # Y|X,Z
              dnorm(t[i], mean = meanXZ, sd = sdXZ)* # X|Z
              pnorm(t[i], mean = meanCXZ[i], sd = sdCXZ, lower.tail = FALSE) # C|X,Z
          }
          return(val)
        }
        
        ### Evaluate
        j_numerator_integrate <- function(jj=1) {
          integrate(function(y) {integral_func_num(t=y,j=jj)}, 
                    lower = -10, upper = 10, rel.tol = 1e-7, abs.tol = 1e-7)$value}
        v.area_num <- Vectorize(j_numerator_integrate)
        numerator = v.area_num(1:3)
        
        # bottom integral
        integral_func_denom <- function(t=1){
          val = rep(NA, length(t))
          for (i in 1:length(t)){
            val[i] = likelihood_int(t[i])*dnorm(t[i], mean = meanXZ, sd = sdXZ) # X|Z
          }
          return(val)
        }
        
        ### Evaluate
        denominator = 
          integrate(integral_func_denom, lower = -10, upper = 10, rel.tol = 1e-7, abs.tol = 1e-7)$value
        
        # Compute Augmented part 
        return(numerator/denominator)
      }
      
      # print(theta)
      ACC_est = D*CCscore() + (1-D/myp)*psi_hat_i()
      return(ACC_est)
    }
  }
  
  # if eta1==TRUE, then wrong augmentation
  if(eta1){
    
    gee_estfun <- function(data, formula){
      
      W <- model.matrix(object = ~ A+Z, data = data)
      X <- model.matrix(object = formula, data = data)
      Y <- model.response(model.frame(formula = formula, data = data))
      D <- data$D
      z <- data$Z
      myp <- data$myp
      
      # need these for integration components
      meanXZ = mymeanXZ(data)
      
      function(theta){
        
        #################
        # 1st Component #
        #################
        CCscore <- function(){
          # useful quantities
          mu <- X %*% theta
          e = Y - mu
          # score equations for betas
          dotbeta = (e/psi^2)%*%X
          return(c(dotbeta))
        }
        
        #######################
        # Augmented Component #
        #######################
        # augmentation component
        psi_hat_i <- function(){
          
          # values needed
          e_star = Y - (W %*% theta)
          a = - 1/(2*sdXZ^2) -theta[2]^2/(2*psi^2) 
          b = (meanXZ/sdXZ^2 -theta[2]*e_star/psi^2)
          c = (-meanXZ^2/(2*sdXZ^2) - e_star^2/(2*psi^2))
          mu_star = -b/(2*a)
          sd2_star = (sdXZ^2*psi^2)/(psi^2 + theta[2]^2*sdXZ^2)
          
          ### beta0
          beta0_star = (theta[2]*mu_star + e_star)/psi^2
          ### beta1
          beta1_star = (-theta[2]*(sd2_star+mu_star^2) + mu_star*(W[2]*theta[2] - e_star) + W[2]*e_star)/psi^2
          ### beta2 
          beta2_star = W[2]*beta0_star
          
          return(c(beta0_star, beta1_star, beta2_star))
        }
        
        # print(theta)
        ACC_est = D*CCscore() + (1-D/myp)*psi_hat_i()
        return(ACC_est)
      }
    }
  }
  
  ######## EST INITIAL PARAMETER ESTIMATES ##########
  mylmer  = lm(model, data = data_yXZ, weight = D)
  g = c(coef(mylmer))
  psi = summary(mylmer)$sigma
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    root_control = setup_root_control(start = c(g)),
    outer_args = list(formula = model),
    deriv_control = setup_deriv_control(method = 'Richardson')
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% matrix()
  data_yXZ = data_yXZ %>% subset(D==1)  
  psi_updated = sum((data_yXZ$y - model.matrix(object = model, data = data_yXZ)%*%beta_estimates)^2/(nrow(data_yXZ))) %>% sqrt
  beta_estimates = rbind(beta_estimates, psi_updated)
  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% matrix()
  return(list(beta_est = t(beta_estimates), se_est = t(se_estimates)))
}

estimate_beta_acc_lambda_close <- function(data_yXZ, model, myweight = "oracle"){
  
  ######################################################
  # Step 1: define the parameters of interest
  sig = 0.5
  mu_x = 0; mu_z = 0; 
  var_x = 1; var_z = 1; var_xz = sig;  
  sdXZ = sqrt(var_x-var_xz^2/var_z)
  
  # define eta1 misspecification #
  mymeanXZ = function(data.){mu_x + (var_xz/var_z)*(data.$Z-mu_z)}
  
  # define weight components #
  if (myweight  == "oracle") data_yXZ$myp = data_yXZ$myp_yz
  if (myweight  == "uniform") data_yXZ$myp = data_yXZ$myp_uniform
  
  ######################################################
  # Step 2: define initial estimates of theta using ipw
  mylmer  = lm(model, data = data_yXZ, weight = D)
  g = c(coef(mylmer))
  psi = summary(mylmer)$sigma
  
  #####################################################
  # Step 3: Calculate Lambda using the values for theta
  calculate_A = function(){
    
    theta = g
    
    ## step a:set data
    calculate_A.b = function(data, part1 = TRUE){
      
      # data = data_mvn() %>% subset(b==2)
      # theta = c(1,2,3,1) +0.1
      W <- model.matrix(object = ~ A + Z, data = data)
      X <- model.matrix(object = model, data = data)
      Y <- model.response(model.frame(formula = model, data = data))
      D <- data$D
      myp <- data$myp
      z <- data$Z
      
      # need these for integration components
      meanXZ = mymeanXZ(data)
      
      #################
      # 1st Component #
      #################
      CCscore <- function(){
        # useful quantities
        mu <- X %*% theta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        return(c(dotbeta))
      }
      
      #######################
      # Augmented Component #
      #######################
      psi_hat_i <- function(){
        
        # values needed
        e_star = Y - (W %*% theta)
        a = - 1/(2*sdXZ^2) -theta[2]^2/(2*psi^2) 
        b = (meanXZ/sdXZ^2 -theta[2]*e_star/psi^2)
        c = (-meanXZ^2/(2*sdXZ^2) - e_star^2/(2*psi^2))
        mu_star = -b/(2*a)
        sd2_star = (sdXZ^2*psi^2)/(psi^2 + theta[2]^2*sdXZ^2)
        
        ### beta0
        beta0_star = (theta[2]*mu_star + e_star)/psi^2
        ### beta1
        beta1_star = (-theta[2]*(sd2_star+mu_star^2) + mu_star*(W[2]*theta[2] - e_star) + W[2]*e_star)/psi^2
        ### beta2 
        beta2_star = W[2]*beta0_star
        
        return(c(beta0_star, beta1_star, beta2_star))
      }
      
      # print(theta)
      if(part1 == TRUE){
        myreturn = matrix(ncol = 1, (D-myp)*D*CCscore() ) %*% psi_hat_i()
        return(myreturn)
      }
      if(part1 == FALSE){
        myreturn = (D-myp)*psi_hat_i()
        myreturn = myreturn%*%t(myreturn)
      }
    }
    
    # estimate first part of A matrix
    part1 = lapply(1:nrow(data_yXZ), function(x) calculate_A.b(data_yXZ[x,]))
    part1 = Reduce("+", part1) 
    
    part2= lapply(1:nrow(data_yXZ), function(x) calculate_A.b(data_yXZ[x,], part1=FALSE))
    part2 = Reduce("+", part2) 
    
    returnA = -part1%*%solve(part2)
    return(returnA)  
  }
  myA = calculate_A()
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    W <- model.matrix(object = ~ A+Z, data = data)
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    D <- data$D
    z <- data$Z
    myp <- data$myp
    
    # need these for integration components
    meanXZ = mymeanXZ(data)
    
    function(theta){
      
      #################
      # 1st Component #
      #################
      CCscore <- function(){
        # useful quantities
        mu <- X %*% theta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        return(c(dotbeta))
      }
      
      #######################
      # Augmented Component #
      #######################
      psi_hat_i <- function(){
        
        # values needed
        e_star = Y - (W %*% theta)
        a = - 1/(2*sdXZ^2) -theta[2]^2/(2*psi^2) 
        b = (meanXZ/sdXZ^2 -theta[2]*e_star/psi^2)
        c = (-meanXZ^2/(2*sdXZ^2) - e_star^2/(2*psi^2))
        mu_star = -b/(2*a)
        sd2_star = (sdXZ^2*psi^2)/(psi^2 + theta[2]^2*sdXZ^2)
        
        ### beta0
        beta0_star = (theta[2]*mu_star + e_star)/psi^2
        ### beta1
        beta1_star = (-theta[2]*(sd2_star+mu_star^2) + mu_star*(W[2]*theta[2] - e_star) + W[2]*e_star)/psi^2
        ### beta2 
        beta2_star = W[2]*beta0_star
        
        return(c(beta0_star, beta1_star, beta2_star))
      }
      
      # print(theta)
      ACC_est = CCscore()*D + (D-myp)*myA%*%psi_hat_i()
      return(ACC_est)
    }
  }
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    root_control = setup_root_control(start = c(g)),
    outer_args = list(formula = model),
    deriv_control = setup_deriv_control(method = 'Richardson')
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% matrix()
  data_yXZ = data_yXZ %>% subset(D==1)  
  psi_updated = sum((data_yXZ$y - model.matrix(object = model, data = data_yXZ)%*%beta_estimates)^2/(nrow(data_yXZ))) %>% sqrt
  beta_estimates = rbind(beta_estimates, psi_updated)
  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% matrix()
  return(list(beta_est = t(beta_estimates), se_est = t(se_estimates)))
}

# macc
estimate_beta_macc <- function(data_yXZ, model, eta1 = FALSE, myweight = "oracle", dep_censoring = FALSE){
  
  ######################################################
  # define parameters values for integration component #
  ######################################################
  sig = 0.5
  mu_x = 0; mu_z = 0; mu_c =0
  var_x = 1; var_z = 1; var_c = 1
  var_xz = sig; var_cz = var_xz; var_xc = var_xz*var_cz 
  sdXZ = sqrt(var_x-var_xz^2/var_z)
  sdCXZ = sqrt(var_c-var_cz^2/var_z)
  
  if (dep_censoring){
    var_xc = 0.6
    sdCXZ = sqrt(var_c - (var_x*var_z - var_xz^2)^(-1)*
                   ((var_xc*var_z-var_cz*var_xz)*(var_xc) +
                      (var_cz*var_x-var_xc*var_xz)*(var_cz)))
  }
  
  ################################
  # define eta1 misspecification #
  ################################
  mymeanXZ = function(data.){mu_x + (var_xz/var_z)*(data.$Z-mu_z)}
  
  if (myweight  == "oracle") data_yXZ$myp = data_yXZ$myp_xz
  if (myweight  == "cox") data_yXZ$myp = data_yXZ$myp_xz_cox
  if (myweight  == "mle") data_yXZ$myp = data_yXZ$myp_xz_mle
  if (myweight  == "uniform") data_yXZ$myp = data_yXZ$myp_uniform
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    W <- model.matrix(object = ~ Z + A, data = data)
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    D <- data$D
    z <- data$Z
    myp <- data$myp
    
    # need these for integration components
    meanXZ = mymeanXZ(data)
    
    function(theta){
      
      #################
      # 1st Component #
      #################
      CCscore <- function(){
        # useful quantities
        mu <- X %*% theta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        return(c(dotbeta))
      }
      
      #######################
      # Augmented Component #
      #######################
      # augmentation component
      psi_hat_i <- function(){
        
        # top integral
        likelihood_int = function(t=1){
          return(dnorm(Y -  cbind(W[,1],W[,2]-t, W[,3]) %*% theta, 0, sd = psi))
        }
        
        # score of integral
        score_int <- function(t=1,j=1){
          # generate empty array of values to save
          # create temporary X matrix
          Wtemp <- cbind(W[,1],W[,2]-t, W[,3])
          mu <- Wtemp %*% theta
          e = Y - mu
          # score equations for theta
          dotbeta = (e/psi^2)%*%Wtemp
          return(c(dotbeta)[j])
        }
        
        #top integral#
        integral_func_num <- function(t=1,j=1){
          val = rep(NA, length(t))
          meanCXZ = mu_c + (var_x*var_z - var_xz^2)^(-1)*
            ((var_xc*var_z-var_cz*var_xz)*(t - mu_x) +
               (var_cz*var_x-var_xc*var_xz)*(z - mu_z))
          for (i in 1:length(t)){
            val[i] = score_int(t[i],j=j)*likelihood_int(t[i])*
              dnorm(t[i], mean = meanXZ, sd = sdXZ)* # X|Z
              (-pnorm(t[i], mean = meanCXZ[i], sd = sdCXZ, lower.tail = TRUE)) # C|X,Z 
          }
          return(val)
        }
        
        ### Evaluate
        j_numerator_integrate <- function(jj=1) {
          integrate(function(y) {integral_func_num(t=y,j=jj)}, 
                    lower = -10, upper = 10, rel.tol = 1e-7, abs.tol = 1e-7)$value}
        v.area_num <- Vectorize(j_numerator_integrate)
        numerator = v.area_num(1:3)
        
        # bottom integral
        integral_func_denom <- function(t=1){
          val = rep(NA, length(t))
          meanCXZ = mu_c + (var_x*var_z - var_xz^2)^(-1)*
            ((var_xc*var_z-var_cz*var_xz)*(t - mu_x) +
               (var_cz*var_x-var_xc*var_xz)*(z - mu_z))
          for (i in 1:length(t)){
            val[i] = likelihood_int(t[i])*
              dnorm(t[i], mean = meanXZ, sd = sdXZ)* # X|Z
              (1/pnorm(t[i], mean = meanCXZ[i], sd = sdCXZ, lower.tail = FALSE)-1) # 1/pi
          }
          return(val)
        }
        
        ### Evaluate
        denominator = 
          integrate(integral_func_denom, lower = -10, upper = 10, rel.tol = 1e-7, abs.tol = 1e-7)$value #-1
        
        # Compute Augmented part 
        return(-numerator/denominator)
      }
      
      # print(theta)
      MACC_est = D*CCscore() + (1-D/myp)*psi_hat_i()
      return(MACC_est)
    }
  }
  
  # if eta1==TRUE, then wrong augmentation
  if(eta1){
    
    gee_estfun <- function(data, formula){
      
      W <- model.matrix(object = ~ Z + A, data = data)
      X <- model.matrix(object = formula, data = data)
      Y <- model.response(model.frame(formula = formula, data = data))
      D <- data$D
      z <- data$Z
      myp <- data$myp
      
      # need these for integration components
      meanXZ = mymeanXZ(data)
      
      function(theta){
        
        #################
        # 1st Component #
        #################
        CCscore <- function(){
          # useful quantities
          mu <- X %*% theta
          e = Y - mu
          # score equations for betas
          dotbeta = (e/psi^2)%*%X
          return(c(dotbeta))
        }
        
        #######################
        # Augmented Component #
        #######################
        # augmentation component
        psi_hat_i <- function(){
          
          # values needed
          e_star = Y - (W %*% theta)
          a = - 1/(2*sdXZ^2) -theta[2]^2/(2*psi^2) 
          b = (meanXZ/sdXZ^2 -theta[2]*e_star/psi^2)
          c = (-meanXZ^2/(2*sdXZ^2) - e_star^2/(2*psi^2))
          mu_star = -b/(2*a)
          sd2_star = (sdXZ^2*psi^2)/(psi^2 + theta[2]^2*sdXZ^2)
          
          ### beta0
          beta0_star = (theta[2]*mu_star + e_star)/psi^2
          ### beta1
          beta1_star = (-theta[2]*(sd2_star+mu_star^2) + mu_star*(W[2]*theta[2] - e_star) + W[2]*e_star)/psi^2
          ### beta2 
          beta2_star = W[2]*beta0_star
          
          return(c(beta0_star, beta1_star, beta2_star))
        }
        
        # print(theta)
        MACC_est = D*CCscore() + (1-D/myp)*psi_hat_i()
        return(MACC_est)
      }
    }
  }
  
  ######## EST INITIAL PARAMETER ESTIMATES ##########
  mylmer  = lm(model, data = data_yXZ, weight = D)
  g = c(coef(mylmer))
  psi = summary(mylmer)$sigma
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    root_control = setup_root_control(start = c(g)),
    outer_args = list(formula = model),
    deriv_control = setup_deriv_control(method = 'Richardson')
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% matrix()
  data_yXZ = data_yXZ %>% subset(D==1)  
  psi_updated = sum((data_yXZ$y - model.matrix(object = model, data = data_yXZ)%*%beta_estimates)^2/(nrow(data_yXZ))) %>% sqrt
  beta_estimates = rbind(beta_estimates, psi_updated)
  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% matrix()
  return(list(beta_est = t(beta_estimates), se_est = t(se_estimates)))
}

estimate_beta_macc_lambda_close <- function(data_yXZ, model, myweight = "oracle"){
  
  ######################################################
  # Step 1: define the parameters of interest
  sig = 0.5
  mu_x = 0; mu_z = 0; 
  var_x = 1; var_z = 1; var_xz = sig;  
  sdXZ = sqrt(var_x-var_xz^2/var_z)
  
  # define eta1 misspecification #
  mymeanXZ = function(data.){mu_x + (var_xz/var_z)*(data.$Z-mu_z)}
  
  # define weight components #
  if (myweight  == "oracle") data_yXZ$myp = data_yXZ$myp_xz
  if (myweight  == "cox") data_yXZ$myp = data_yXZ$myp_xz_cox
  if (myweight  == "mle") data_yXZ$myp = data_yXZ$myp_xz_mle
  if (myweight  == "uniform") data_yXZ$myp = data_yXZ$myp_uniform
  
  ######################################################
  # Step 2: define initial estimates of theta using ipw
  mylmer  = lm(model, data = data_yXZ, weight = D)
  g = c(coef(mylmer))
  psi = summary(mylmer)$sigma
  
  #####################################################
  # Step 3: Calculate Lambda using the values for theta
  calculate_A = function(){
    
    theta = g
    
    ## step a:set data
    calculate_A.b = function(data, part1 = TRUE){
      
      # data = data_mvn() %>% subset(b==2)
      # theta = c(1,2,3,1) +0.1
      W <- model.matrix(object = ~ A+Z, data = data)
      X <- model.matrix(object = model, data = data)
      Y <- model.response(model.frame(formula = model, data = data))
      D <- data$D
      myp <- data$myp
      z <- data$Z
      
      # need these for integration components
      meanXZ = mymeanXZ(data)
      
      #################
      # 1st Component #
      #################
      CCscore <- function(){
        # useful quantities
        mu <- X %*% theta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        return(c(dotbeta))
      }
      
      #######################
      # Augmented Component #
      #######################
      psi_hat_i <- function(){
        
        # values needed
        e_star = Y - (W %*% theta)
        a = - 1/(2*sdXZ^2) -theta[2]^2/(2*psi^2) 
        b = (meanXZ/sdXZ^2 -theta[2]*e_star/psi^2)
        c = (-meanXZ^2/(2*sdXZ^2) - e_star^2/(2*psi^2))
        mu_star = -b/(2*a)
        sd2_star = (sdXZ^2*psi^2)/(psi^2 + theta[2]^2*sdXZ^2)
        
        ### beta0
        beta0_star = (theta[2]*mu_star + e_star)/psi^2
        ### beta1
        beta1_star = (-theta[2]*(sd2_star+mu_star^2) + mu_star*(W[2]*theta[2] - e_star) + W[2]*e_star)/psi^2
        ### beta2 
        beta2_star = W[2]*beta0_star
        
        return(c(beta0_star, beta1_star, beta2_star))
      }
      
      # print(theta)
      if(part1 == TRUE){
        myreturn = matrix(ncol = 1, (1-D/myp)*D*CCscore() ) %*% psi_hat_i()
        return(myreturn)
      }
      if(part1 == FALSE){
        myreturn = psi_hat_i()
        myreturn = (1-D/myp)^2*myreturn%*%t(myreturn)
      }
    }
    
    # estimate first part of A matrix
    part1 = lapply(1:nrow(data_yXZ), function(x) calculate_A.b(data_yXZ[x,]))
    part1 = Reduce("+", part1) 
    
    part2= lapply(1:nrow(data_yXZ), function(x) calculate_A.b(data_yXZ[x,], part1=FALSE))
    part2 = Reduce("+", part2) 
    
    returnA = -part1%*%solve(part2)
    return(returnA)  
  }
  myA = calculate_A()
  #print("done with A"); print(myA)
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    W <- model.matrix(object = ~ A+Z, data = data)
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    D <- data$D
    z <- data$Z
    myp <- data$myp
    
    # need these for integration components
    meanXZ = mymeanXZ(data)
    
    function(theta){
      
      #################
      # 1st Component #
      #################
      CCscore <- function(){
        # useful quantities
        mu <- X %*% theta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        return(c(dotbeta))
      }
      
      #######################
      # Augmented Component #
      #######################
      psi_hat_i <- function(){
        
        # values needed
        e_star = Y - (W %*% theta)
        a = - 1/(2*sdXZ^2) -theta[2]^2/(2*psi^2) 
        b = (meanXZ/sdXZ^2 -theta[2]*e_star/psi^2)
        c = (-meanXZ^2/(2*sdXZ^2) - e_star^2/(2*psi^2))
        mu_star = -b/(2*a)
        sd2_star = (sdXZ^2*psi^2)/(psi^2 + theta[2]^2*sdXZ^2)
        
        ### beta0
        beta0_star = (theta[2]*mu_star + e_star)/psi^2
        ### beta1
        beta1_star = (-theta[2]*(sd2_star+mu_star^2) + mu_star*(W[2]*theta[2] - e_star) + W[2]*e_star)/psi^2
        ### beta2 
        beta2_star = W[2]*beta0_star
        
        return(c(beta0_star, beta1_star, beta2_star))
      }
      
      # print(theta)
      MACC_est = CCscore()*D + (1 - D/myp)* myA%*%psi_hat_i()
      return(MACC_est)
    }
  }
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    root_control = setup_root_control(start = c(g)),
    outer_args = list(formula = model),
    deriv_control = setup_deriv_control(method = 'Richardson')
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% matrix()
  data_yXZ = data_yXZ %>% subset(D==1)  
  psi_updated = sum((data_yXZ$y - model.matrix(object = model, data = data_yXZ)%*%beta_estimates)^2/(nrow(data_yXZ))) %>% sqrt
  beta_estimates = rbind(beta_estimates, psi_updated)
  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% matrix()
  return(list(beta_est = t(beta_estimates), se_est = t(se_estimates)))
}

# aipw
estimate_beta_aipw <- function(data_yXZ, model, eta1 = FALSE, myweight = "oracle", dep_censoring=FALSE){
  
  ######################################################
  # define parameters values for integration component #
  ######################################################
  sig = 0.5
  mu_x = 0; mu_z = 0; mu_c = 0
  var_x = 1; var_z = 1; var_c = 1
  var_xz = sig; var_cz = var_xz; var_xc = var_xz*var_cz 
  sdXZ = sqrt(var_x-var_xz^2/var_z)
  sdCXZ = sqrt(var_c-var_cz^2/var_z)
  
  if (dep_censoring){
    var_xc = 0.6
    sdCXZ = sqrt(var_c - (var_x*var_z - var_xz^2)^(-1)*
                   ((var_xc*var_z-var_cz*var_xz)*(var_xc) +
                      (var_cz*var_x-var_xc*var_xz)*(var_cz)))
  }
  
  ################################
  # define eta1 misspecification #
  ################################
  mymeanXZ = function(data.){mu_x + (var_xz/var_z)*(data.$Z-mu_z)}
  # if(eta1){mymeanXZ = function(data.){-2}}
  
  if (myweight  == "oracle") data_yXZ$myp = data_yXZ$myp_xz
  if (myweight  == "cox") data_yXZ$myp = data_yXZ$myp_xz_cox
  if (myweight  == "mle") data_yXZ$myp = data_yXZ$myp_xz_mle
  if (myweight  == "uniform") data_yXZ$myp = data_yXZ$myp_uniform
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    W <- model.matrix(object = ~ A+Z, data = data)
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    D <- data$D
    z <- data$Z
    myp <- data$myp
    
    # need these for integration components
    meanXZ = mymeanXZ(data)
    
    # meanCXZ // needs to be estimated for each case of X
    
    function(theta){
      
      #################
      # 1st Component #
      #################
      CCscore <- function(){
        # useful quantities
        mu <- X %*% theta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        return(c(dotbeta))
      }
      
      #######################
      # Augmented Component #
      #######################
      # augmentation component
      psi_hat_i <- function(){
        
        # top integral
        likelihood_int = function(t=1){
          return(dnorm(Y -  cbind(W[,1],W[,2]-t, W[,3]) %*% theta, 0, sd = psi))
        }
        
        # score of integral
        score_int <- function(t=1,j=1){
          # generate empty array of values to save
          # create temporary X matrix
          Wtemp <- cbind(W[,1],W[,2]-t, W[,3])
          mu <- Wtemp %*% theta
          e = Y - mu
          # score equations for theta
          dotbeta = (e/psi^2)%*%Wtemp
          return(c(dotbeta)[j])
        }
        
        #top integral#
        ## switch tt instead of t
        integral_func_num <- function(t=1,j=1){
          val = rep(NA, length(t))
          meanCXZ = mu_c + (var_x*var_z - var_xz^2)^(-1)*
            ((var_xc*var_z-var_cz*var_xz)*(t - mu_x) +
               (var_cz*var_x-var_xc*var_xz)*(z - mu_z))
          for (i in 1:length(t)){
            val[i] = score_int(t[i],j=j)*likelihood_int(t[i])*
              dnorm(t[i], mean = meanXZ, sd = sdXZ)* # X|Z
              (1-1/pnorm(t[i], mean = meanCXZ[i], sd = sdCXZ, lower.tail = FALSE)) # (1-1/pi)
          }
          return(val)
        }
        
        ### Evaluate
        j_numerator_integrate <- function(jj=1) {
          integrate(function(y) {integral_func_num(t=y,j=jj)}, 
                    lower = -10, upper = 10, rel.tol = 1e-7, abs.tol = 1e-7)$value}
        v.area_num <- Vectorize(j_numerator_integrate)
        numerator = v.area_num(1:3)
        
        # bottom integral
        integral_func_denom <- function(t=1){
          val = rep(NA, length(t))
          meanCXZ = mu_c + (var_x*var_z - var_xz^2)^(-1)*
            ((var_xc*var_z-var_cz*var_xz)*(t - mu_x) +
               (var_cz*var_x-var_xc*var_xz)*(z - mu_z))
          for (i in 1:length(t)){
            val[i] = likelihood_int(t[i])*
              dnorm(t[i], mean = meanXZ, sd = sdXZ)* # X|Z
              (1-1/pnorm(t[i], mean = meanCXZ[i], sd = sdCXZ, lower.tail = FALSE)) # (1-1/pi)
          }
          return(val)
        }
        
        ### Evaluate
        denominator = integrate(integral_func_denom, 
                                lower = -10, upper = 10, rel.tol = 1e-7, abs.tol = 1e-7)$value #+ 1
        
        # Compute Augmented part 
        return(numerator/denominator)
      }
      

      return(CCscore()*D/myp + (1 - D/myp)*psi_hat_i())
    }
  }
  
  if(eta1){
    
    gee_estfun <- function(data, formula){
      
      W <- model.matrix(object = ~ A + Z, data = data)
      X <- model.matrix(object = formula, data = data)
      Y <- model.response(model.frame(formula = formula, data = data))
      D <- data$D
      z <- data$Z
      myp <- data$myp
      
      # need these for integration components
      meanXZ = mymeanXZ(data)
      
      function(theta){
        
        #################
        # 1st Component #
        #################
        CCscore <- function(){
          # useful quantities
          mu <- X %*% theta
          e = Y - mu
          # score equations for betas
          dotbeta = (e/psi^2)%*%X
          return(c(dotbeta))
        }
        
        #######################
        # Augmented Component #
        #######################
        psi_hat_i <- function(){
          
          # values needed
          e_star = Y - (W %*% theta)
          a = - 1/(2*sdXZ^2) -theta[2]^2/(2*psi^2) 
          b = (meanXZ/sdXZ^2 -theta[2]*e_star/psi^2)
          c = (-meanXZ^2/(2*sdXZ^2) - e_star^2/(2*psi^2))
          mu_star = -b/(2*a)
          sd2_star = (sdXZ^2*psi^2)/(psi^2 + theta[2]^2*sdXZ^2)
          
          ### beta0
          beta0_star = (theta[2]*mu_star + e_star)/psi^2
          ### beta1
          beta1_star = (-theta[2]*(sd2_star+mu_star^2) + mu_star*(W[2]*theta[2] - e_star) + W[2]*e_star)/psi^2
          ### beta2 
          beta2_star = W[2]*beta0_star
          
          return(c(beta0_star, beta1_star, beta2_star))
        }
        
        # print(theta)
        AIPW_est = D*CCscore()/myp + (1-D/myp)*psi_hat_i()
        return(AIPW_est)
      }
    }
  }
  
  ######## EST INITIAL PARAMETER ESTIMATES ##########
  mylmer  = lm(model, data = data_yXZ, weights = D/myp)
  g = c(coef(mylmer))
  psi = summary(mylmer)$sigma
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    root_control = setup_root_control(start = c(g)),
    outer_args = list(formula = model),
    deriv_control = setup_deriv_control(method = 'Richardson')
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% matrix()
  data_yXZ = data_yXZ %>% subset(D==1)  
  psi_updated = sum((data_yXZ$y - model.matrix(object = model, data = data_yXZ)%*%beta_estimates)^2/(nrow(data_yXZ))) %>% sqrt
  beta_estimates = rbind(beta_estimates, psi_updated)
  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% matrix()
  return(list(beta_est = t(beta_estimates), se_est = t(se_estimates)))
}

estimate_beta_aipw_lambda_close <- function(data_yXZ, model, myweight = "oracle"){
  
  ######################################################
  # Step 1: define the parameters of interest
  sig = 0.5
  mu_x = 0; mu_z = 0; 
  var_x = 1; var_z = 1; var_xz = sig;  
  sdXZ = sqrt(var_x-var_xz^2/var_z)
  
  # define eta1 misspecification #
  mymeanXZ = function(data.){mu_x + (var_xz/var_z)*(data.$Z-mu_z)}
  
  #######################################################
  # define weight components #
  #######################################################
  if (myweight  == "oracle") data_yXZ$myp = data_yXZ$myp_xz
  if (myweight  == "cox") data_yXZ$myp = data_yXZ$myp_xz_cox
  if (myweight  == "mle") data_yXZ$myp = data_yXZ$myp_xz_mle
  if (myweight  == "uniform") data_yXZ$myp = data_yXZ$myp_uniform
  
  ######################################################
  # Step 2: define initial estimates of theta using ipw
  mylmer  = lm(model, data = data_yXZ, weight = D/myp)
  g = c(coef(mylmer))
  psi = summary(mylmer)$sigma
  
  #####################################################
  # Step 3: Calculate Lambda using the values for theta
  calculate_A = function(){
    
    theta = g
    
    ## step a:set data
    calculate_A.b = function(data, part1 = TRUE){
      
      W <- model.matrix(object = ~ A+Z, data = data)
      X <- model.matrix(object = model, data = data)
      Y <- model.response(model.frame(formula = model, data = data))
      D <- data$D
      myp <- data$myp
      z <- data$Z
      
      # need these for integration components
      meanXZ = mymeanXZ(data)
      
      #################
      # 1st Component #
      #################
      CCscore <- function(){
        # useful quantities
        mu <- X %*% theta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        return(c(dotbeta))
      }
      
      #######################
      # Augmented Component #
      #######################
      # augmentation component
      
      psi_hat_i <- function(){
        
        # values needed
        e_star = Y - (W %*% theta)
        a = - 1/(2*sdXZ^2) -theta[2]^2/(2*psi^2) 
        b = (meanXZ/sdXZ^2 -theta[2]*e_star/psi^2)
        c = (-meanXZ^2/(2*sdXZ^2) - e_star^2/(2*psi^2))
        mu_star = -b/(2*a)
        sd2_star = (sdXZ^2*psi^2)/(psi^2 + theta[2]^2*sdXZ^2)
        
        ### beta0
        beta0_star = (theta[2]*mu_star + e_star)/psi^2
        ### beta1
        beta1_star = (-theta[2]*(sd2_star+mu_star^2) + mu_star*(W[2]*theta[2] - e_star) + W[2]*e_star)/psi^2
        ### beta2 
        beta2_star = W[2]*beta0_star
        
        return(c(beta0_star, beta1_star, beta2_star))
      }
     
      # print(theta)
      if(part1 == TRUE){
        myreturn = psi_hat_i()
        myreturn = matrix(ncol = 1, (1-D/myp)*(D/myp)*CCscore()) %*% psi_hat_i()
        return(myreturn)
      }
      if(part1 == FALSE){
        myreturn = psi_hat_i()
        myreturn = (1-D/myp)^2*myreturn%*%t(myreturn)
      }
    }
    
    # estimate first part of A matrix
    part1 = lapply(1:nrow(data_yXZ), function(x) calculate_A.b(data_yXZ[x,], part1=TRUE))
    part1 = Reduce("+", part1) 
    
    part2 = lapply(1:nrow(data_yXZ), function(x) calculate_A.b(data_yXZ[x,], part1=FALSE))
    part2 = Reduce("+", part2) 
    
    returnA = -part1%*%solve(part2)
    return(returnA)  
  }
  myA = calculate_A()
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    W <- model.matrix(object = ~ A + Z, data = data)
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    D <- data$D
    z <- data$Z
    myp <- data$myp
    
    # need these for integration components
    meanXZ = mymeanXZ(data)
    
    # meanCXZ // needs to be estimated for each case of X
    
    function(theta){
      
      #################
      # 1st Component #
      #################
      CCscore <- function(){
        # useful quantities
        mu <- X %*% theta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        return(c(dotbeta))
      }
      
      #######################
      # Augmented Component #
      #######################
      # augmentation component
      psi_hat_i <- function(){
        
        # values needed
        e_star = Y - (W %*% theta)
        a = - 1/(2*sdXZ^2) -theta[2]^2/(2*psi^2) 
        b = (meanXZ/sdXZ^2 -theta[2]*e_star/psi^2)
        c = (-meanXZ^2/(2*sdXZ^2) - e_star^2/(2*psi^2))
        mu_star = -b/(2*a)
        sd2_star = (sdXZ^2*psi^2)/(psi^2 + theta[2]^2*sdXZ^2)
        
        ### beta0
        beta0_star = (theta[2]*mu_star + e_star)/psi^2
        ### beta1
        beta1_star = (-theta[2]*(sd2_star+mu_star^2) + mu_star*(W[2]*theta[2] - e_star) + W[2]*e_star)/psi^2
        ### beta2 
        beta2_star = W[2]*beta0_star
        
        return(c(beta0_star, beta1_star, beta2_star))
      }
      
      # print(theta)
      AIPW_est = CCscore()*D/myp + (1 - D/myp)* myA %*% psi_hat_i()
      return(AIPW_est)
    }
  }
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    root_control = setup_root_control(start = c(g)),
    outer_args = list(formula = model),
    deriv_control = setup_deriv_control(method = 'Richardson')
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% matrix()
  data_yXZ = data_yXZ %>% subset(D==1)  
  psi_updated = sum((data_yXZ$y - model.matrix(object = model, data = data_yXZ)%*%beta_estimates)^2/(nrow(data_yXZ))) %>% sqrt
  beta_estimates = rbind(beta_estimates, psi_updated)
  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% matrix()
  return(list(beta_est = t(beta_estimates), se_est = t(se_estimates)))
}

# mle

estimate_beta_likelihood <- function(data_yXZ, model, eta1 = FALSE, dep_censoring=FALSE){
  
  ######################################################
  # define parameters values for integration component #
  ######################################################
  sig = 0.5
  mu_x = 0; mu_z = 0; mu_c =0
  var_x = 1; var_z = 1; var_c = 1
  var_xz = sig; var_cz = var_xz; var_xc = var_xz*var_cz 
  sdXZ = sqrt(var_x-var_xz^2/var_z)
  sdCXZ = sqrt(var_c-var_cz^2/var_z)
  
  if (dep_censoring){
    var_xc = 0.6
    sdCXZ = sqrt(var_c - (var_x*var_z - var_xz^2)^(-1)*
                   ((var_xc*var_z-var_cz*var_xz)*(var_xc) +
                      (var_cz*var_x-var_xc*var_xz)*(var_cz)))
  }
  
  ################################
  # define eta1 misspecification #
  ################################
  data_yXZ$meanXZ = mu_x + (var_xz/var_z)*(data_yXZ$Z-mu_z)
  if(eta1){data_yXZ$meanXZ = -2}
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    W <- model.matrix(object = ~ A + Z, data = data)
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    D <- data$D
    w <- data$W
    z <- data$Z
    
    # need these for integration components
    meanXZ =  data$meanXZ
    
    function(theta){
      
      # set theta
      beta = theta[1:3]
      psi = theta[4]
      
      # first component of score
      CCscore <- function(){
        # useful quantities
        mu <- X %*% beta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        # score equations for sigmas
        dotsigma = -(1/psi) + e^2/psi^3
        return(c(dotbeta,dotsigma))
      }
      if(D ==1) return(CCscore())
      
      # second component of score
      e_myscore <- function(){
        
        # top integral
        likelihood_int = function(t=1){
          return(dnorm(Y -  cbind(W[,1],W[,2]-t, W[,3]) %*% beta, 0, sd = psi))
        }
        
        # score of integral
        score_int <- function(t=1,j=1){
          # generate empty array of values to save
          # create temporary X matrix
          Wtemp <- cbind(W[,1],W[,2]-t, W[,3])
          mu <- Wtemp %*% beta
          e = Y - mu
          # score equations for betas
          dotbeta = (e/psi^2)%*%Wtemp
          # score equations for sigmas
          dotsigma = -(1/psi) + e^2/psi^3
          return(c(dotbeta,dotsigma)[j])
        }
        
        #top integral#
        integral_func_psi <- function(t=1,j=1){
          val = rep(NA, length(t))
          meanCXZ = mu_c + (var_x*var_z - var_xz^2)^(-1)*
            ((var_xc*var_z-var_cz*var_xz)*(t - mu_x) +
               (var_cz*var_x-var_xc*var_xz)*(z - mu_z))
          for (i in 1:length(t)){
            val[i] = score_int(t[i],j=j)*
              likelihood_int(t[i])*
              dnorm(t[i], mean =  meanXZ, sd = sdXZ)* # X|Z
              dnorm(w, mean =  meanCXZ[i], sd = sdCXZ)
          }
          return(val)
        }
        
        ## evaluate
        area_num <- function(j) integrate(integral_func_psi,j=j,
                                          lower = w,
                                          upper = 10)$value
        v.area_num <- Vectorize(area_num)
        numerator = v.area_num(1:4)
        
        # bottom integral
        integral_func_denom <- function(t=1){
          val = rep(NA, length(t))
          meanCXZ = mu_c + (var_x*var_z - var_xz^2)^(-1)*
            ((var_xc*var_z-var_cz*var_xz)*(t - mu_x) +
               (var_cz*var_x-var_xc*var_xz)*(z - mu_z))
          for (i in 1:length(t)){
            val[i] = likelihood_int(t[i])*
              dnorm(t[i], mean =  meanXZ, sd = sdXZ)* # X|Z
              dnorm(w, mean =  meanCXZ[i], sd = sdCXZ)
          }
          return(val)
        }
        # evaluate
        denominator = integrate(function(x) integral_func_denom(t=x),
                                lower = w,
                                upper = 10)$value
        
        # evaluate and return
        return(numerator/denominator)
      }
      return(e_myscore())
      
    
      
    }
  }
  
  ######## EST INITIAL PARAMETER ESTIMATES ##########
  mylmer  = lm(model, data = data_yXZ %>% subset(D==1))
  g = c(coef(mylmer), summary(mylmer)$sigma)
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    root_control = setup_root_control(start = c(g)),
    outer_args = list(formula = model)
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% as.data.frame()
  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% as.data.frame()
  return(list(beta_est = beta_estimates, se_est = se_estimates[1:3]))
}

#################################
# Estimated nuisance parameters #
#################################

#ipw 
var_beta_ipw <- function(data_yXZ, mytheta, myalpha){
  
  # data_yXZ = data_yXZ
  # mytheta = unlist(est30$beta_est)
  # myalpha = myalpha1
  
  mybeta = mytheta[1:(length(mytheta)-1)]
  mypsi = mytheta[length(mytheta)]
  myxi = c(mybeta,myalpha)
  
  #########################################################
  # alpha helper functions
  pi_xz_func.b = function(data., myalpha.){
    pnorm(data.$W, 
          mean = myalpha.[3:4]%*%c(1, data.$Z),
          sd = myalpha.[6], lower.tail = FALSE)
  }
  
  alpha.logLik.b = function(myalpha., data.){
    
    # density function
    cxz_likelihood = function(w,delta){
      
      if(delta==0){
        mymeanXZ = myalpha.[1:2]%*%c(1, data.$Z)
        mymeanCZ = myalpha.[3:4]%*%c(1, data.$Z)
        return(pnorm(w, mymeanXZ, myalpha.[5], lower.tail = FALSE)*
                 dnorm(w,mymeanCZ, myalpha.[6]))
      }
      
      if(delta==1){
        mymeanXZ = myalpha.[1:2]%*%c(1, data.$Z)
        mymeanCZ = myalpha.[3:4]%*%c(1, data.$Z)
        return(pnorm(w, mymeanCZ, myalpha.[6], lower.tail = FALSE)*
                 dnorm(w,mymeanXZ, myalpha.[5]))
      }
      
    }
    
    # integrate the density function above
    myreturn = cxz_likelihood(data.$W, data.$D) 
    return(myreturn)
  }
  alpha.jacobian.b = function(data){
    myderiv = numDeriv::jacobian(function(x) 
      logLik.b(myalpha. = x, data. = data), myalpha, method = "Richardson")
    return(matrix(ncol=1, myderiv))
  }
  alpha.hessian.b = function(data){
    myderiv = numDeriv::hessian(function(x) 
      logLik.b(myalpha. = x, data. = data), myalpha, method = "Richardson")
    return(myderiv)
  }
  
  #########################################################
  # define the stacked score equation
  ipwscore.b <- function(myxi., data.b){
    
    # myxi. = myxi; data.b = dat[1,]
    
    mybeta. = myxi.[1:length(mybeta)]
    myalpha. = myxi.[(length(mybeta)+1):length(myxi.)]
    
    ######################
    ## IPW
    X <- model.matrix(object =  y~AW +Z, data = data.b)
    Y <- model.response(model.frame(formula =  y~AW+Z, data = data.b))
    myp = pi_xz_func.b(data.=data.b, myalpha. = myalpha.)
    D <- data.b$D
    mu <- X %*% mybeta.
    e = Y - mu
    # score equations for betas
    dotbeta = (e/mypsi^2)%*%X
    dotbeta = D*c(dotbeta)/myp
    
    ########################
    ## alpha 
    dotalpha = numDeriv::jacobian(function(x) 
      alpha.logLik.b(myalpha. = x, data. = data.b), myalpha., method = "Richardson")
    
    return(c(dotbeta,dotalpha))
  }
  
  # A matrix
  calculate.A = function(){
    part1 = parallel::mclapply(1:nrow(data_yXZ), function(i) 
      dotalpha = numDeriv::jacobian(function(x) 
        ipwscore.b(myxi. = x, data.b = data_yXZ[i,]), myxi, method = "Richardson")
    )
    part1 = Reduce("+", part1) 
  }
  myA = calculate.A()
  myA.inv = solve(myA)
  
  # B matrix
  calculate.B = function(){
    part1 = parallel::mclapply(1:nrow(data_yXZ), function(i){
      dotalpha = ipwscore.b(myxi. = myxi, data.b = data_yXZ[i,])
      return(dotalpha %*% t(dotalpha))} 
    )
    part1 = Reduce("+", part1)
    return(part1)
  }
  myB = calculate.B()
  sand.var = myA.inv%*%myB%*%t(myA.inv)
  # return my sandwich variance estimator 
  
  return(list(se_updated = sqrt(diag(sand.var))[1:3]))  
}

# macc
var_beta_macc_lambda <- function(data_yXZ, mytheta, myalpha){
  
  # dat$myp = dat$myp_xz
  # data_yXZ = dat
  # mytheta = unlist(est60$beta_est)
  # myalpha = myalpha1
  
  ######################################################
  # Step 1: Calculate Lambda using the values for theta
  g = mytheta[1:(length(mytheta)-1)]
  psi = mytheta[length(mytheta)]
  
  sdXZ = var(data_yXZ$meanXZ_hat)^0.5
  meanXZ = mean(data_yXZ$meanXZ_hat)
  
  data_yXZ$myp = data_yXZ$myp_xz
  #####################################################
  # Step 2: Define helper functions
  pi_xz_func.b = function(data., myalpha.){
    pnorm(data.$W, 
          mean = myalpha.[3:4]%*%c(1, data.$Z),
          sd = myalpha.[6], lower.tail = FALSE)
  }
  
  jacobian_pi.b = function(data){
    myderiv = numDeriv::jacobian(function(x) 
      pi_xz_func.b(myalpha. = x, data. = data), myalpha, method = "Richardson")
    return(matrix(ncol=1, myderiv))
  } 
  
  
  logLik.b = function(myalpha., data.){
    
    # density function
    cxz_likelihood = function(w,delta){
      
      if(delta==0){
        mymeanXZ = myalpha.[1:2]%*%c(1, data.$Z)
        mymeanCZ = myalpha.[3:4]%*%c(1, data.$Z)
        return(pnorm(w, mymeanXZ, myalpha.[5], lower.tail = FALSE)*
                 dnorm(w,mymeanCZ, myalpha.[6]))
      }
      
      if(delta==1){
        mymeanXZ = myalpha.[1:2]%*%c(1, data.$Z)
        mymeanCZ = myalpha.[3:4]%*%c(1, data.$Z)
        return(pnorm(w, mymeanCZ, myalpha.[6], lower.tail = FALSE)*
                 dnorm(w,mymeanXZ, myalpha.[5]))
      }
      
    }
    
    # integrate the density function above
    myreturn = cxz_likelihood(data.$W, data.$D) 
    return(myreturn)
  }
  
  jacobian.b = function(data){
    myderiv = numDeriv::jacobian(function(x) 
      logLik.b(myalpha. = x, data. = data), myalpha, method = "Richardson")
    return(matrix(ncol=1, myderiv))
  }
  
  hessian.b = function(data){
    myderiv = numDeriv::hessian(function(x) 
      logLik.b(myalpha. = x, data. = data), myalpha, method = "Richardson")
    return(myderiv)
  }
  
  #####################################################
  ## 3. The A_alpha matrix 
  calculate_A_alpha = function(){
    # calculate the hessian matrix and return
    part1 = parallel::mclapply(1:nrow(data_yXZ), function(x) hessian.b(data_yXZ[x,]))
    part1 = Reduce("+", part1) /nrow(data_yXZ)
    return(part1)
  }
  my_A_alpha_inv = -solve(calculate_A_alpha())
  
  ## 2. A_alpha_theta  
  calculate_A_alpha_theta = function(){
    
    theta = g
    
    # augmentation component
    phi_yz.b = function(data.){
      
      W <- model.matrix(object = ~ A + Z, data = data.)
      X <- model.matrix(object = y ~ AW + Z, data = data.)
      Y <- model.response(model.frame(formula = y ~ AW + Z, data = data.))
      D <- data.$D
      myp <- data.$myp
      
      #######################
      # Augmented Component #
      #######################
      psi_hat_i <- function(){
        
        # values needed
        e_star = Y - (W %*% theta)
        a = - 1/(2*sdXZ^2) -theta[2]^2/(2*psi^2) 
        b = (meanXZ/sdXZ^2 -theta[2]*e_star/psi^2)
        c = (-meanXZ^2/(2*sdXZ^2) - e_star^2/(2*psi^2))
        mu_star = -b/(2*a)
        sd2_star = (sdXZ^2*psi^2)/(psi^2 + theta[2]^2*sdXZ^2)
        
        # # denominator
        # mydenom = exp(c-(b^2/(4*a)))/sqrt(2*pi*(psi^2+theta[3]^2*sdXZ^2))
        
        ### beta0
        beta0_star = (theta[2]*mu_star + e_star)/psi^2
        ### beta1
        beta1_star = (-theta[2]*(sd2_star+mu_star^2) + mu_star*(W[2]*theta[2] - e_star) + W[2]*e_star)/psi^2
        ### beta2
        beta2_star = W[3]*beta0_star
        # ### beta3
        # beta3_star = W[4]*beta0_star
        # ### beta4
        # beta4_star = W[5]*beta0_star
        
        return(c(beta0_star, beta1_star, beta2_star))
      }
      
      return(psi_hat_i())
    }
    
    # calculate the value for the ith observations
    A_alpha_theta.b = function(mydata){
      D = mydata$D
      myp = mydata$myp
      matrix(ncol=1, D*phi_yz.b(mydata)/(myp^2)) %*% t(jacobian_pi.b(mydata))
    }
    
    part1 = lapply(1:nrow(data_yXZ), function(x) A_alpha_theta.b(data_yXZ[x,]))
    part1 = Reduce("+", part1) /nrow(data_yXZ)
    return(part1)
  }
  my_A_alpha_theta = calculate_A_alpha_theta()
  
  myconstant = my_A_alpha_theta%*%my_A_alpha_inv
  # print(myconstant)
  #####################################################
  # Step 4: Calculate Lambda using the values for theta
  
  # score function for alpha
  calculate_A = function(){
    
    theta = g
    
    ## calculate the bth component
    calculate_A.b = function(data, part1 = TRUE){
      
      W <- model.matrix(object = ~ A + Z, data = data)
      X <- model.matrix(object = ~ AW + Z, data = data)
      Y <- model.response(model.frame(formula = y ~ AW + Z, data = data))
      D <- data$D
      myp <- data$myp
      
      #################
      # 1st Component #
      #################
      CCscore <- function(){
        # useful quantities
        mu <- X %*% theta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        return(c(dotbeta))
      }
      
      #######################
      # Augmented Component #
      #######################
      psi_hat_i <- function(){
        
        # values needed
        e_star = Y - (W %*% theta)
        a = - 1/(2*sdXZ^2) -theta[2]^2/(2*psi^2) 
        b = (meanXZ/sdXZ^2 -theta[2]*e_star/psi^2)
        c = (-meanXZ^2/(2*sdXZ^2) - e_star^2/(2*psi^2))
        mu_star = -b/(2*a)
        sd2_star = (sdXZ^2*psi^2)/(psi^2 + theta[2]^2*sdXZ^2)
        
        
        ### beta0
        beta0_star = (theta[2]*mu_star + e_star)/psi^2
        ### beta1
        beta1_star = (-theta[2]*(sd2_star+mu_star^2) + mu_star*(W[2]*theta[2] - e_star) + W[2]*e_star)/psi^2
        ### beta2
        beta2_star = W[3]*beta0_star
        # ### beta3
        # beta3_star = W[4]*beta0_star
        # ### beta4
        # beta4_star = W[5]*beta0_star
        
        return(c(beta0_star, beta1_star, beta2_star))
      }
      
      # print(theta)
      if(part1 == TRUE){
        myreturn = matrix(ncol = 1,D*CCscore()) %*% 
          t( (1-D/myp)*matrix(ncol=1,psi_hat_i()) + myconstant%*%jacobian.b(data))
        return(myreturn)
      }
      if(part1 == FALSE){
        myreturn = (1-D/myp)*matrix(ncol=1,psi_hat_i()) + myconstant%*%jacobian.b(data)
        myreturn = myreturn%*%t(myreturn)
      }
    }
    
    # estimate first part of A matrix
    part1 = parallel::mclapply(1:nrow(data_yXZ), function(x) calculate_A.b(data_yXZ[x,]))
    part1 = Reduce("+", part1) 
    
    part2 = parallel::mclapply(1:nrow(data_yXZ), function(x) calculate_A.b(data_yXZ[x,], part1=FALSE))
    part2 = Reduce("+", part2) 
    
    returnA = -part1%*%solve(part2)
    return(returnA)  
  }
  myA = calculate_A()
  print("A is done"); print(myA)
  
  #####################################################
  # Step 5: Define estimating equation
  acc.lambda.b <- function(theta, data.b){
    
    W <- model.matrix(object = ~ A + Z, data = data.b)
    X <- model.matrix(object = y ~ AW + Z, data = data.b)
    Y <- model.response(model.frame(formula = y ~ AW + Z, data = data.b))
    D <- data.b$D
    myp <- data.b$myp
    
    #################
    # 1st Component #
    #################
    CCscore <- function(){
      # useful quantities
      mu <- X %*% theta
      e = Y - mu
      # score equations for betas
      dotbeta = (e/psi^2)%*%X
      return(c(dotbeta))
    }
    
    #######################
    # Augmented Component #
    #######################
    psi_hat_i <- function(){
      
      # values needed
      e_star = Y - (W %*% theta)
      a = - 1/(2*sdXZ^2) -theta[2]^2/(2*psi^2) 
      b = (meanXZ/sdXZ^2 -theta[2]*e_star/psi^2)
      c = (-meanXZ^2/(2*sdXZ^2) - e_star^2/(2*psi^2))
      mu_star = -b/(2*a)
      sd2_star = (sdXZ^2*psi^2)/(psi^2 + theta[2]^2*sdXZ^2)
      
      ### beta0
      beta0_star = (theta[2]*mu_star + e_star)/psi^2
      ### beta1
      beta1_star = (-theta[2]*(sd2_star+mu_star^2) + mu_star*(W[2]*theta[2] - e_star) + W[2]*e_star)/psi^2
      ### beta2
      beta2_star = W[3]*beta0_star
      # ### beta3
      # beta3_star = W[4]*beta0_star
      # ### beta4
      # beta4_star = W[5]*beta0_star
      
      return(c(beta0_star, beta1_star, beta2_star))
    }
    
    # print(theta)
    ACC_est = CCscore()*D + myA %*% ((1 - D/myp)*psi_hat_i() + myconstant%*%jacobian.b(data.b))
    # print(ACC_est)
    return(ACC_est)
    
  }
  
  #####################################################
  # Step 6: Calculate the A and the B matrices
  
  calculate.B = function(){
    part1 = parallel::mclapply(1:nrow(data_yXZ), function(i){
      dotalpha = acc.lambda.b(theta = g, data.b = data_yXZ[i,])
      return(dotalpha %*% t(dotalpha))} 
    )
    part1 = Reduce("+", part1)
    return(part1)
  }
  myB = calculate.B()
  
  # A matrix
  calculate.A = function(){
    part1 = parallel::mclapply(1:nrow(data_yXZ), function(i) 
      dotalpha = numDeriv::jacobian(function(x) 
        acc.lambda.b(theta = x, data.b = data_yXZ[i,]), g, method = "Richardson")
    )
    part1 = Reduce("+", part1) #/nrow(data_yXZ)
  }
  myA = calculate.A()
  myA.inv = solve(myA)
  
  # calculate sandwich variance estimate 
  print("bread"); print(myA.inv)
  print("meat"); print(myB)
  
  sand.var = myA.inv%*%myB%*%t(myA.inv)
  sand.var = sqrt(diag(sand.var))
  
  # return values
  return(list(se_updated = sand.var))
}

# aipw
var_beta_aipw_lambda <- function(data_yXZ, mytheta, myalpha){
  
  ######################################################
  # Step 1: Calculate Lambda using the values for theta
  g = mytheta[1:(length(mytheta)-1)]
  psi = mytheta[length(mytheta)]
  
  sdXZ = var(data_yXZ$meanXZ_hat)^0.5
  meanXZ = mean(data_yXZ$meanXZ_hat)
  
  data_yXZ$myp = data_yXZ$myp_xz
  
  #####################################################
  # Step 2: Define helper functions
  pi_xz_func.b = function(data., myalpha.){
    pnorm(data.$W, 
          mean = myalpha.[3:4]%*%c(1, data.$Z),
          sd = myalpha.[6], lower.tail = FALSE)
  }
  
  jacobian_pi.b = function(data){
    myderiv = numDeriv::jacobian(function(x) 
      pi_xz_func.b(myalpha. = x, data. = data), myalpha, method = "Richardson")
    return(matrix(ncol=1, myderiv))
  } 
  
  logLik.b = function(myalpha., data.){
    
    # density function
    cxz_likelihood = function(w,delta){
      
      if(delta==0){
        mymeanXZ = myalpha.[1:2]%*%c(1, data.$Z)
        mymeanCZ = myalpha.[3:4]%*%c(1, data.$Z)
        return(pnorm(w, mymeanXZ, myalpha.[5], lower.tail = FALSE)*
                 dnorm(w,mymeanCZ, myalpha.[6]))
      }
      
      if(delta==1){
        mymeanXZ = myalpha.[1:2]%*%c(1, data.$Z)
        mymeanCZ = myalpha.[3:4]%*%c(1, data.$Z)
        return(pnorm(w, mymeanCZ, myalpha.[6], lower.tail = FALSE)*
                 dnorm(w,mymeanXZ, myalpha.[5]))
      }
      
    }
    
    # integrate the density function above
    myreturn = cxz_likelihood(data.$W, data.$D) 
    return(myreturn)
  }
  
  jacobian.b = function(data){
    myderiv = numDeriv::jacobian(function(x) 
      logLik.b(myalpha. = x, data. = data), myalpha, method = "Richardson")
    return(matrix(ncol=1, myderiv))
  }
  
  hessian.b = function(data){
    myderiv = numDeriv::hessian(function(x) 
      logLik.b(myalpha. = x, data. = data), myalpha, method = "Richardson")
    return(myderiv)
  }
  
  ## 1. The A_alpha matrix 
  calculate_A_alpha = function(){
    # calculate the hessian matrix and return
    part1= parallel::mclapply(1:nrow(data_yXZ), function(x) hessian.b(data_yXZ[x,]))
    part1 = Reduce("+", part1)/nrow(data_yXZ)
    return(part1)
  }
  my_A_alpha_inv = -solve(calculate_A_alpha())
  
  ## 2. The A_ipw of theta
  calculate_A_ipw = function(){
    
    theta = g
    
    A_ipw.b = function(data) {
      
      X <- model.matrix(object = ~ AW+Z, data = data)
      Y <- model.response(model.frame(formula = y~ AW+Z, data = data))
      D <- data$D
      
      # calculate the IPW score equation
      IPWscore <- function(myalpha.){
        # useful quantities
        mu = X %*% theta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        myp = pi_xz_func.b(data.=data, myalpha. = myalpha.)
        myreturn = matrix(ncol=1, D*c(dotbeta)/myp)
        return(myreturn)
      }
      
      # calculate the partial derivative
      myderiv = numDeriv::jacobian(function(x) 
        IPWscore(myalpha. = x), myalpha, method = "Richardson")
      return(myderiv)
    }
    
    # calculate the hessian matrix and return
    part1= parallel::mclapply(1:nrow(data_yXZ), function(x) A_ipw.b(data_yXZ[x,]))
    part1 = Reduce("+", part1)/nrow(data_yXZ)
    return(part1)
  }
  my_A_ipw = calculate_A_ipw()
  
  ## 3. A_alpha_theta  
  calculate_A_alpha_theta = function(){
    
    theta = g
    
    # augmentation component
    phi_yz.b = function(data.){
      
      W <- model.matrix(object = ~ A + Z, data = data.)
      X <- model.matrix(object = y ~ AW + Z, data = data.)
      Y <- model.response(model.frame(formula = y ~ AW + Z, data = data.))
      D <- data.$D
      myp <- data.$myp
      z <- data.$Z
      
      #######################
      # Augmented Component #
      #######################
      psi_hat_i <- function(){
        
        # values needed
        e_star = Y - (W %*% theta)
        a = - 1/(2*sdXZ^2) -theta[2]^2/(2*psi^2) 
        b = (meanXZ/sdXZ^2 -theta[2]*e_star/psi^2)
        c = (-meanXZ^2/(2*sdXZ^2) - e_star^2/(2*psi^2))
        mu_star = -b/(2*a)
        sd2_star = (sdXZ^2*psi^2)/(psi^2 + theta[2]^2*sdXZ^2)
        
        ### beta0
        beta0_star = (theta[2]*mu_star + e_star)/psi^2
        ### beta1
        beta1_star = (-theta[2]*(sd2_star+mu_star^2) + mu_star*(W[2]*theta[2] - e_star) + W[2]*e_star)/psi^2
        ### beta2
        beta2_star = W[3]*beta0_star
        
        return(c(beta0_star, beta1_star, beta2_star))
      }
      
      return(psi_hat_i())
    }
    
    # calculate the value for the ith observations
    A_alpha_theta.b = function(mydata){
      D = mydata$D
      myp = mydata$myp
      matrix(ncol=1, D*phi_yz.b(mydata)/(myp^2)) %*% t(jacobian_pi.b(mydata))
    }
    
    part1= parallel::mclapply(1:nrow(data_yXZ), function(x) A_alpha_theta.b(data_yXZ[x,]))
    part1 = Reduce("+", part1) /nrow(data_yXZ)
    return(part1)
  }
  my_A_alpha_theta = calculate_A_alpha_theta()
  
  #####################################################
  # Step 4: Calculate Lambda using the values for theta
  
  # score function for alpha
  calculate_A = function(){
    
    theta = g
    
    ## calculate the bth component
    calculate_A.b = function(data, part1 = TRUE){
      
      # data = data_mvn() %>% subset(b==2)
      # theta = c(1,2,3) 
      W <- model.matrix(object = ~ A + Z, data = data)
      X <- model.matrix(object = y ~ AW + Z, data = data)
      Y <- model.response(model.frame(formula = y ~ AW + Z, data = data))
      D <- data$D
      myp <- data$myp
      
      #################
      # 1st Component #
      #################
      CCscore <- function(){
        # useful quantities
        mu <- X %*% theta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        return(c(dotbeta))
      }
      
      #######################
      # Augmented Component #
      #######################
      psi_hat_i <- function(){
        
        # values needed
        e_star = Y - (W %*% theta)
        a = - 1/(2*sdXZ^2) -theta[2]^2/(2*psi^2) 
        b = (meanXZ/sdXZ^2 -theta[2]*e_star/psi^2)
        c = (-meanXZ^2/(2*sdXZ^2) - e_star^2/(2*psi^2))
        mu_star = -b/(2*a)
        sd2_star = (sdXZ^2*psi^2)/(psi^2 + theta[2]^2*sdXZ^2)
        
        ### beta0
        beta0_star = (theta[2]*mu_star + e_star)/psi^2
        ### beta1
        beta1_star = (-theta[2]*(sd2_star+mu_star^2) + mu_star*(W[2]*theta[2] - e_star) + W[2]*e_star)/psi^2
        ### beta2
        beta2_star = W[3]*beta0_star
        
        return(c(beta0_star, beta1_star, beta2_star))
      }
      
      # print(theta)
      if(part1 == TRUE){
        myreturn = matrix(ncol=1, (1-D/myp)*psi_hat_i()) + my_A_alpha_theta%*%my_A_alpha_inv%*%jacobian.b(data)
        myreturn = (matrix(ncol = 1,D*CCscore()/myp) + my_A_ipw%*%my_A_alpha_inv%*%jacobian.b(data)) %*% t(myreturn) 
        return(myreturn)
      }
      if(part1 == FALSE){
        myreturn = matrix(ncol=1, (1-D/myp)*psi_hat_i()) + my_A_alpha_theta%*%my_A_alpha_inv%*%jacobian.b(data)
        myreturn = myreturn%*%t(myreturn)
      }
    }
    
    # estimate first part of A matrix
    part1 = parallel::mclapply(1:nrow(data_yXZ), function(x) calculate_A.b(data_yXZ[x,]))
    part1 = Reduce("+", part1)/nrow(data_yXZ) 
    
    part2= lapply(1:nrow(data_yXZ), function(x) calculate_A.b(data_yXZ[x,], part1=FALSE))
    part2 = Reduce("+", part2)/nrow(data_yXZ) 
    
    returnA = -part1%*%solve(part2)
    return(returnA)  
  }
  myA = calculate_A()
  print("A is done"); print(myA)
  toc()
  
  #####################################################
  # Step 6: Calculate the A and the B matrices
  aipw.lambda.b <- function(theta, data.b){
    
    W <- model.matrix(object = ~ A + Z, data = data.b)
    X <- model.matrix(object = y ~ AW + Z, data = data.b)
    Y <- model.response(model.frame(formula = y ~ AW + Z, data = data.b))
    D <- data.b$D
    myp <- data.b$myp
    
    #################
    # 1st Component #
    #################
    CCscore <- function(){
      # useful quantities
      mu <- X %*% theta
      e = Y - mu
      # score equations for betas
      dotbeta = (e/psi^2)%*%X
      return(c(dotbeta))
    }
    
    #######################
    # Augmented Component #
    #######################
    
    psi_hat_i <- function(){
      
      # values needed
      e_star = Y - (W %*% theta)
      a = - 1/(2*sdXZ^2) -theta[2]^2/(2*psi^2) 
      b = (meanXZ/sdXZ^2 -theta[2]*e_star/psi^2)
      c = (-meanXZ^2/(2*sdXZ^2) - e_star^2/(2*psi^2))
      mu_star = -b/(2*a)
      sd2_star = (sdXZ^2*psi^2)/(psi^2 + theta[2]^2*sdXZ^2)
      
      ### beta0
      beta0_star = (theta[2]*mu_star + e_star)/psi^2
      ### beta1
      beta1_star = (-theta[2]*(sd2_star+mu_star^2) + mu_star*(W[2]*theta[2] - e_star) + W[2]*e_star)/psi^2
      ### beta2
      beta2_star = W[3]*beta0_star
      
      return(c(beta0_star, beta1_star, beta2_star))
    }
    
    # print(theta)
    AIPW_est = CCscore()*D/myp + (1 - D/myp)* myA %*% psi_hat_i() + 
      (my_A_ipw +  myA %*% my_A_alpha_theta)%*%(my_A_alpha_inv%*%jacobian.b(data.b))
    return(AIPW_est)
    
  }
  
  #####################################################
  # Step 6: Calculate the A and the B matrices
  
  calculate.B = function(){
    part1 = parallel::mclapply(1:nrow(data_yXZ), function(i){
      dotalpha = aipw.lambda.b(theta = g, data.b = data_yXZ[i,])
      return(dotalpha %*% t(dotalpha))} 
    )
    part1 = Reduce("+", part1)
    return(part1)
  }
  myB = calculate.B()
  
  # A matrix
  calculate.A = function(){
    part1 = parallel::mclapply(1:nrow(data_yXZ), function(i) 
      dotalpha = numDeriv::jacobian(function(x) 
        aipw.lambda.b(theta = x, data.b = data_yXZ[i,]), g, method = "Richardson")
    )
    part1 = Reduce("+", part1)
  }
  myA = calculate.A()
  myA.inv = solve(myA)
  
  # calculate sandwich variance estimate 
  sand.var = myA.inv%*%myB%*%t(myA.inv)
  sand.var = sqrt(diag(sand.var))
  # return values
  return(list(se_updated = sand.var))
}

# mle
estimate_beta_likelihood_hat <- function(data_yXZ, model, myalpha){
  
  ######################################################
  # define parameters values for integration component #
  ######################################################
  sdXZ = myalpha[length(myalpha)-1]
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    W <- model.matrix(object = ~ A+Z, data = data)
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    D <- data$D
    w <- data$W
    z <- data$Z
    b <- data$b
    meanXZ = data$meanXZ_hat
    
    function(theta){
      
      # set theta
      beta = theta[1:3]
      psi = theta[4]
      
      # first component of score
      CCscore <- function(){
        # useful quantities
        mu <- X %*% beta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        dotsigma = -(1/psi) + e^2/psi^3
        return(c(dotbeta, dotsigma))
      }
      
      # second component of score
      e_myscore <- function(){
        
        # top integral
        likelihood_int = function(t=1){
          return(Y-dnorm(cbind(W[,1],W[,2]-t, W[,3]) %*% beta, 0, sd = psi))
        }
        
        # score of integral
        score_int <- function(t=1,j=1){
          # generate empty array of values to save
          # create temporary X matrix
          Wtemp <- cbind(W[,1],W[,2]-t, W[,3])
          e = Y - Wtemp %*% beta
          # score equations for betas
          dotbeta = (e/psi^2)%*%Wtemp
          dotsigma = -(1/psi) + e^2/psi^3
          return(c(dotbeta,dotsigma)[j])
        }
        
        #top integral#
        integral_func_psi <- function(t=1,j=1){
          val = rep(NA, length(t))
          for (i in 1:length(t)){
            val[i] = score_int(t[i],j=j)*likelihood_int(t[i])*
              dnorm(t[i], mean =  meanXZ, sd = sdXZ) 
          }
          return(val)
        }
        
        ## evaluate
        area_num <- function(j) integrate(integral_func_psi,j=j,
                                          lower = w,
                                          upper = Inf)$value
        v.area_num <- Vectorize(area_num)
        numerator = v.area_num(1:4)
        
        # bottom integral
        integral_func_denom <- function(t=1){
          val = rep(NA, length(t))
          for (i in 1:length(t)){
            val[i] = likelihood_int(t[i])*
              dnorm(t[i], mean =  meanXZ, sd = sdXZ) 
          }
          return(val)
        }
        # evaluate
        denominator = integrate(function(x) integral_func_denom(t=x),
                                lower = w,
                                upper = Inf)$value
        
        # evaluate and return
        return(numerator/denominator)
      }
      
      # return imputation
      # print(theta); print(b)
      return(D*CCscore() + (1-D)*e_myscore())
    }
  }
  
  ######## EST INITIAL PARAMETER ESTIMATES ##########
  mylmer  = lm(model, data = data_yXZ, weights = D)
  g = c(coef(mylmer), summary(mylmer)$sigma)
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    root_control = setup_root_control(start = c(g)),
    outer_args = list(formula = model)
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% as.data.frame()
  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% as.data.frame()
  return(list(beta_est = beta_estimates, se_est = se_estimates[1:3]))
}

var_beta_likelihood <- function(data_yXZ, mytheta, myalpha){
  
  ######################################################
  # Step 1: Calculate Lambda using the values for theta 
  mytheta = unlist(mytheta)
  sdXZ = myalpha[length(myalpha)-1]
  
  #####################################################
  # Step 2: Define helper functions
  logLik.b = function(myalpha., data.){
    
    # density function
    cxz_likelihood = function(w,delta){
      
      if(delta==0){
        mymeanXZ = myalpha.[1:2]%*%c(1, data.$Z)
        mymeanCZ = myalpha.[3:4]%*%c(1, data.$Z)
        return(pnorm(w, mymeanXZ, myalpha.[5], lower.tail = FALSE)*
                 dnorm(w,mymeanCZ, myalpha.[6]))
      }
      
      if(delta==1){
        mymeanXZ = myalpha.[1:2]%*%c(1, data.$Z)
        mymeanCZ = myalpha.[3:4]%*%c(1, data.$Z)
        return(pnorm(w, mymeanCZ, myalpha.[6], lower.tail = FALSE)*
                 dnorm(w,mymeanXZ, myalpha.[5]))
      }
      
    }
    
    # integrate the density function above
    myreturn = cxz_likelihood(data.$W, data.$D) 
    return(myreturn)
  }
  
  jacobian.b = function(data){
    myderiv = numDeriv::jacobian(function(x) 
      logLik.b(myalpha. = x, data. = data), myalpha, method = "Richardson")
    return(matrix(ncol=1, myderiv))
  }
  
  hessian.b = function(data){
    myderiv = numDeriv::hessian(function(x) 
      logLik.b(myalpha. = x, data. = data), myalpha, method = "Richardson")
    return(myderiv)
  }
  
  #####################################################
  ## 3. The A_alpha matrix 
  calculate_A_alpha = function(){
    # calculate the hessian matrix and return
    part1 = parallel::mclapply(1:nrow(data_yXZ), function(x) hessian.b(data_yXZ[x,]))
    part1 = Reduce("+", part1) /nrow(data_yXZ)
    return(part1)
  }
  my_A_alpha_inv = -solve(calculate_A_alpha())
  
  ## 2. A_alpha_theta  
  calculate_A_alpha_theta = function(){
    
    beta = mytheta[1:3]
    psi = mytheta[4]
    # augmentation component
    phi_yz.b = function(myalpha., data.){
      
      W <- model.matrix(object = y ~ A + Z, data = data.)
      X <- model.matrix(object = y ~ AW + Z, data = data.)
      Y <- model.response(model.frame(formula = y ~ AW + Z, data = data.))
      D <- data.$D
      w <- data.$W
      meanXZ = myalpha.[1:2]%*%c(1, data.$Z)
      sdXZ = myalpha.[5]
      
      e_myscore_i <- function(){
        
        # top integral
        likelihood_int = function(t=1){
          return(Y-dnorm(cbind(W[,1],W[,2]-t, W[,3]) %*% beta, 0, sd = psi))
        }
        
        # score of integral
        score_int <- function(t=1,j=1){
          # generate empty array of values to save
          # create temporary X matrix
          Wtemp <- cbind(W[,1],W[,2]-t, W[,3])
          e = Y - Wtemp %*% beta
          # score equations for betas
          dotbeta = (e/psi^2)%*%Wtemp
          # score equations for sigmas
          dotsigma = -(1/psi) + e^2/psi^3
          return(c(dotbeta,dotsigma)[j])
        }
        
        #top integral#
        integral_func_psi <- function(t=1,j=1){
          val = rep(NA, length(t))
          for (i in 1:length(t)){
            val[i] = score_int(t[i],j=j)*likelihood_int(t[i])*
              dnorm(t[i], mean =  meanXZ, sd = sdXZ) 
          }
          return(val)
        }
        
        ## evaluate
        area_num <- function(j) integrate(integral_func_psi,j=j,
                                          lower = w,
                                          upper = Inf)$value
        v.area_num <- Vectorize(area_num)
        numerator = v.area_num(1:4)
        
        # bottom integral
        integral_func_denom <- function(t=1){
          val = rep(NA, length(t))
          for (i in 1:length(t)){
            val[i] = likelihood_int(t[i])*
              dnorm(t[i], mean =  meanXZ, sd = sdXZ) 
          }
          return(val)
        }
        # evaluate
        denominator = integrate(function(x) integral_func_denom(t=x),
                                lower = w,
                                upper = Inf)$value
        
        # evaluate and return
        return(numerator/denominator)
      }
      
      return(e_myscore_i())
    }
    
    part1 = parallel::mclapply(1:nrow(data_yXZ), function(i) 
      dotalpha = numDeriv::jacobian(function(x) 
        (1-data_yXZ$D[i])*phi_yz.b(myalpha. = x, data. = data_yXZ[i,]), myalpha, method = "Richardson"))
    
    part1 = Reduce("+", part1) /nrow(data_yXZ)
    return(part1)
  }
  my_A_alpha_theta = calculate_A_alpha_theta()
  myconstant = my_A_alpha_theta%*%my_A_alpha_inv
  
  #####################################################
  # Step 5: Define estimating equation
  likelihood.lambda.b <- function(theta, data.b){
    
    W <- model.matrix(object = ~ A + Z, data = data.b)
    X <- model.matrix(object = y ~ AW + Z, data = data.b)
    Y <- model.response(model.frame(formula = y ~ AW + Z, data = data.b))
    D <- data.b$D
    w <- data.b$W
    meanXZ = data.b$meanXZ_hat
    
    # set theta
    beta = theta[1:3]
    psi = theta[4]
    
    # first component of score
    CCscore <- function(){
      # useful quantities
      mu <- X %*% beta
      e = Y - mu
      # score equations for betas
      dotbeta = (e/psi^2)%*%X
      # score equations for sigmas
      dotsigma = -(1/psi) + e^2/psi^3
      return(c(dotbeta, dotsigma))
    }
    
    # second component of score
    e_myscore <- function(){
      
      # top integral
      likelihood_int = function(t=1){
        return(Y-dnorm(cbind(W[,1],W[,2]-t, W[,3]) %*% beta, 0, sd = psi))
      }
      
      # score of integral
      score_int <- function(t=1,j=1){
        # generate empty array of values to save
        # create temporary X matrix
        Wtemp <- cbind(W[,1],W[,2]-t, W[,3])
        e = Y - Wtemp %*% beta
        # score equations for betas
        dotbeta = (e/psi^2)%*%Wtemp
        # score equations for sigmas
        dotsigma = -(1/psi) + e^2/psi^3
        return(c(dotbeta, dotsigma)[j])
      }
      
      #top integral#
      integral_func_psi <- function(t=1,j=1){
        val = rep(NA, length(t))
        for (i in 1:length(t)){
          val[i] = score_int(t[i],j=j)*likelihood_int(t[i])*
            dnorm(t[i], mean =  meanXZ, sd = sdXZ) 
        }
        return(val)
      }
      
      ## evaluate
      area_num <- function(j) integrate(integral_func_psi,j=j,
                                        lower = w,
                                        upper = Inf)$value
      v.area_num <- Vectorize(area_num)
      numerator = v.area_num(1:4)
      
      # bottom integral
      integral_func_denom <- function(t=1){
        val = rep(NA, length(t))
        for (i in 1:length(t)){
          val[i] = likelihood_int(t[i])*
            dnorm(t[i], mean =  meanXZ, sd = sdXZ) 
        }
        return(val)
      }
      # evaluate
      denominator = integrate(function(x) integral_func_denom(t=x),
                              lower = w,
                              upper = Inf)$value
      
      # evaluate and return
      return(numerator/denominator)
    }
    
    # return imputation
    return(D*CCscore() + (1-D)*e_myscore())
    
  }
  
  #####################################################
  # Step 6: Calculate the A and the B matrices
  
  calculate.B = function(){
    part1 = parallel::mclapply(1:nrow(data_yXZ), function(i){
      dotalpha = likelihood.lambda.b(theta = mytheta, data.b = data_yXZ[i,])
      return(dotalpha %*% t(dotalpha))} 
    )
    part1 = Reduce("+", part1)
    return(part1)
  }
  myB = calculate.B()
  
  # A matrix
  calculate.A = function(){
    part1 = parallel::mclapply(1:nrow(data_yXZ), function(i) 
      dotalpha = numDeriv::jacobian(function(x) 
        likelihood.lambda.b(theta = x, data.b = data_yXZ[i,]), mytheta, method = "Richardson")
    )
    part1 = Reduce("+", part1)
  }
  myA = calculate.A()
  myA.inv = solve(myA)
  
  # calculate sandwich variance estimate 
  sand.var = myA.inv%*%myB%*%t(myA.inv)
  sand.var = sqrt(diag(sand.var))[1:3]
  # return values
  return(list(se_updated = sand.var))
}

################################
# simulation helper functions # 
###############################

simcalc = function(sim_i=1, dep_censoring. = FALSE){
  set.seed(sim_i)
  dat <- data_mvn(nSubjects = n, dep_censoring = dep_censoring.)
  dat$myp_uniform = runif(n, min=0.1, max = 0.9)
  
  # Oracle 
  est00 = estimate_beta(data_yXZ = dat, model = y ~ AX + Z) # oracle
  
  # Naive
  est10 = estimate_beta(data_yXZ = dat, model = y ~ AW + Z) # naive
  
  # CC
  est20 = estimate_beta_cc(data_yXZ = dat, model = y ~ AW + Z) # CC
  
  # IPW
  est30 = estimate_beta_ipw(data_yXZ = dat, model = y ~ AW+Z, myweight = "oracle") # oracle weight pi(w,z)
  est31 = estimate_beta_ipw(data_yXZ = dat, model = y ~ AW+Z, myweight = "yz") # oracle weight pi(y,z)
  est32 = estimate_beta_ipw(data_yXZ = dat, model = y ~ AW+Z, myweight = "uniform") # random weights
  
  # MLE
  est40 = estimate_beta_likelihood(dat, y ~ AW+Z, eta1= FALSE, dep_censoring = dep_censoring.) # correct specification of f(X,Z)
  est41 = estimate_beta_likelihood(dat, y ~ AW+Z, eta1= TRUE, dep_censoring = dep_censoring.) # incorrect specification of f(X|Z)
  
  # ACC
  est50 = estimate_beta_acc(dat, y ~ AW+Z, eta1= FALSE, myweight = "oracle", dep_censoring = dep_censoring.)
  est51 = estimate_beta_acc(dat, y ~ AW+Z, eta1= TRUE, myweight = "oracle", dep_censoring = dep_censoring.)
  est52 = estimate_beta_acc(dat, y ~ AW+Z, eta1= FALSE, myweight = "uniform", dep_censoring = dep_censoring.)
  est53 = estimate_beta_acc(dat, y ~ AW+Z, eta1= TRUE, myweight = "uniform", dep_censoring = dep_censoring.)
  est54 = estimate_beta_acc_lambda_close(dat, y ~ AW+Z, myweight = "oracle")
  est55 = estimate_beta_acc_lambda_close(dat, y ~ AW+Z, myweight = "uniform")
  
  # MACC
  est60 = estimate_beta_macc(dat, y ~ AW+Z, eta1= FALSE, myweight = "oracle", dep_censoring = dep_censoring.)
  est61 = estimate_beta_macc(dat, y ~ AW+Z, eta1= TRUE, myweight = "oracle", dep_censoring = dep_censoring.)
  est62 = estimate_beta_macc(dat, y ~ AW+Z, eta1= FALSE, myweight = "uniform", dep_censoring = dep_censoring.)
  est63 = estimate_beta_macc(dat, y ~ AW+Z, eta1= TRUE, myweight = "uniform", dep_censoring = dep_censoring.)
  est64 = estimate_beta_macc_lambda_close(dat, y ~ AW+Z, myweight = "oracle")
  est65 = estimate_beta_macc_lambda_close(dat, y ~ AW+Z, myweight = "uniform")
  
  # AIPW
  est70 = estimate_beta_aipw(dat, y ~ AW+Z, eta1= FALSE, myweight = "oracle", dep_censoring = dep_censoring.)
  est71 = estimate_beta_aipw(dat, y ~ AW+Z, eta1= TRUE, myweight = "oracle", dep_censoring = dep_censoring.)
  est72 = estimate_beta_aipw(dat, y ~ AW+Z, eta1= FALSE, myweight = "uniform", dep_censoring = dep_censoring.)
  est73 = estimate_beta_aipw(dat, y ~ AW+Z, eta1= TRUE, myweight = "uniform", dep_censoring = dep_censoring.)
  est74 = estimate_beta_aipw_lambda_close(dat, y ~ AW+Z, myweight = "oracle")
  est75 = estimate_beta_aipw_lambda_close(dat, y ~ AW+Z, myweight = "uniform")
  
  # save results
  myreturn = data.frame(cbind(c("oracle","naive", "cc", 
                                "ipw-wz", "ipw-yz", "ipw-w", 
                                "acc", "acc-psi", "acc-w", "acc-psi-w", "acc-lambda", "acc-lambda-w", 
                                "macc", "macc-psi", "macc-w", "macc-psi-w", "macc-lambda", "macc-lambda-w", 
                                "aipw", "aipw-psi", "aipw-w", "aipw-psi-w", "aipw-lambda", "aipw-lambda-w"),
                              # beta estimates + sigma
                              rbind(est00$beta_est,est10$beta_est,est20$beta_est,
                                    est30$beta_est, est31$beta_est, est32$beta_est, 
                                    est50$beta_est, est51$beta_est, est52$beta_est, est53$beta_est, est54$beta_est, est55$beta_est,
                                    est60$beta_est, est61$beta_est, est62$beta_est, est63$beta_est, est64$beta_est, est65$beta_est,
                                    est70$beta_est, est71$beta_est, est62$beta_est, est73$beta_est, est74$beta_est, est75$beta_est),
                              # standard error estimates
                              rbind(est00$se_est,est10$se_est,est20$se_est,
                                    est30$se_est,est31$se_est,est32$se_est, 
                                    est50$se_est, est51$se_est, est52$se_est, est53$se_est, est54$se_est, est55$se_est,
                                    est60$se_est, est61$se_est, est62$se_est, est63$se_est, est64$se_est, est65$se_est,
                                    est70$se_est, est71$se_est, est72$se_est, est73$se_est, est74$se_est, est75$se_est)
  ))
  colnames(myreturn) = c("Method",
                         "beta0.est", "beta1.est", "beta2.est","sigma.est",
                         "beta0.se", "beta1.se", "beta2.se")
  
  
  myreturn_mle = data.frame(cbind(c("mle","mle-fxz"),
                                  rbind(est40$beta_est, est41$beta_est),
                                  rbind(est40$se_est, est41$se_est)))
  colnames(myreturn_mle) = c("Method",
                             "beta0.est", "beta1.est", "beta2.est","sigma.est",
                             "beta0.se", "beta1.se", "beta2.se") 
  
  # print results 
  myreturn = rbind(myreturn, myreturn_mle)
  return(myreturn)
}