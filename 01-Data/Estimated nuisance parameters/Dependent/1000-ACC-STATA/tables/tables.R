rm(list = ls())
library(tidyverse)
library(xtable)
library(readxl)
`%!in%` = Negate(`%in%`)

#myfolder = "Dependent/ACC"
#myn = "1000-ACC-Test"
mysig = 95
zval  = 1-(1-mysig/100)/2

# Estimates for theta
#list_est <- list.files(path = paste0("../",myfolder,"/",myn), pattern = "*csv", full.names = T)

## Misspecified model
list_est <- list.files(path = paste0("../"), pattern = "augacc_z_mysim*", full.names = T)
myresults_aipw <-  list_est %>% map_df(function(x) read_excel(x, col_names = FALSE))
myresults_est <- myresults_aipw[4*(1:length(list_est))-3, ]
myresults_se <- myresults_aipw[-(4*(1:length(list_est))-3), ]
myresults_se_temp <- matrix(nrow=1, ncol=3)
for (i in 1:length(list_est)){
    temp = myresults_se[(i*3-2):(i*3),]
    myresults_se_temp = rbind(myresults_se_temp, 
        cbind(as.numeric(temp[1,1]), as.numeric(temp[2,2]), as.numeric(temp[3,3])))
}
myresults_se = myresults_se_temp[-1,]^0.5
theta = c(0,0.2,0.2,0.95)
myresults_aipw_z = cbind(myresults_est, myresults_se)
myresults_aipw_z$Method = "acc-z" 

## Correct model
list_est <- list.files(path = paste0("../"), pattern = "augacc_x z_mysim*", full.names = T)
myresults_aipw <-  list_est %>% map_df(function(x) read_excel(x, col_names = FALSE))
myresults_est <- myresults_aipw[4*(1:length(list_est))-3, ]
myresults_se <- myresults_aipw[-(4*(1:length(list_est))-3), ]
myresults_se_temp <- matrix(nrow=1, ncol=3)
for (i in 1:length(list_est)){
    temp = myresults_se[(i*3-2):(i*3),]
    myresults_se_temp = rbind(myresults_se_temp, 
        cbind(as.numeric(temp[1,1]), as.numeric(temp[2,2]), as.numeric(temp[3,3])))
}
myresults_se = myresults_se_temp[-1,]^0.5
theta = c(0,0.2,0.2,0.95)
myresults_aipw_yz = cbind(myresults_est, myresults_se)
myresults_aipw_yz$Method = "acc-yz" 

## combine both models
myresults_aipw = rbind(myresults_aipw_yz, myresults_aipw_z)
colnames(myresults_aipw) = c("beta1.x", "beta2.x", "beta0.x", "beta1.y", "beta2.y", "beta0.y", "Method")

# based on m-estimate function
myresults_aipw$beta0.coverage <- ((myresults_aipw$beta0.x + qnorm(zval)*myresults_aipw$beta0.y) >= theta[1]) &
  ((myresults_aipw$beta0.x - qnorm(zval)*myresults_aipw$beta0.y) <= theta[1])
myresults_aipw$beta1.coverage <- ((myresults_aipw$beta1.x + qnorm(zval)*myresults_aipw$beta1.y) >= theta[2]) &
  ((myresults_aipw$beta1.x - qnorm(zval)*myresults_aipw$beta1.y) <= theta[2])
myresults_aipw$beta2.coverage <- ((myresults_aipw$beta2.x + qnorm(zval)*myresults_aipw$beta2.y) >= theta[3]) &
  ((myresults_aipw$beta2.x - qnorm(zval)*myresults_aipw$beta2.y) <= theta[3])

myresults_est = myresults_aipw %>%
  # select only variables of interest
  dplyr::select(Method, beta0.x, beta1.x, beta2.x) %>%
  # calculate means
  group_by(Method) %>%
  summarise_all(function(x) round(mean(x),3))

myresults_est$beta0.bias = (myresults_est$beta0.x-theta[1])
myresults_est$beta1.bias = (myresults_est$beta1.x-theta[2])/theta[2]
myresults_est$beta2.bias = (myresults_est$beta2.x-theta[3])/theta[3]

myresults_se = myresults_aipw %>%
  # select only variables of interest
  dplyr::select(Method, beta0.y, beta1.y, beta2.y) %>%
  # calculate means
  group_by(Method) %>%
  summarise_all(function(x) round(100*mean(x),3))
  
myresults_sehat = myresults_aipw %>%
  # select only variables of interest
  dplyr::select(Method, beta0.x, beta1.x, beta2.x) %>%
  # calculate the means
  group_by(Method) %>%
  summarise_all(function(x) round(100*var(x)^0.5,3))
colnames(myresults_sehat) = c("Method","beta0.z","beta1.z","beta2.z")

myresults_coverage_mest = myresults_aipw %>%
  # calculate coverage
  mutate(beta0.coverage = ifelse(beta0.coverage, 1, 0),
         beta1.coverage = ifelse(beta1.coverage, 1, 0),
         beta2.coverage = ifelse(beta2.coverage, 1, 0)) %>%
  # select only variables of interest
  dplyr::select(Method, beta0.coverage, beta1.coverage, beta2.coverage) %>%
  # calculate means
  group_by(Method) %>%
  summarise_all(function(x) round(100*mean(x),2))



myresults_out = cbind(myresults_est,myresults_se[,-1], myresults_sehat[,-1], myresults_coverage_mest[,-1])
outputtable = cbind(myresults_out[,1],
                    myresults_out[ , grepl( "beta0" , names(myresults_out))],
                    myresults_out[ , grepl( "beta1" , names(myresults_out))],
                    myresults_out[ , grepl( "beta2" , names(myresults_out))])

print(outputtable)
