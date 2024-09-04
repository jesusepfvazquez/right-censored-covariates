rm(list = ls())
library(tidyverse)
library(xtable)
library(latex2exp)
library(ggplot2)
`%!in%` = Negate(`%in%`)


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
myresults_aipw_z$Method = "acc-lambda-eta1-w" 

## Correct model
list_est <- list.files(path = paste0("../"), pattern = "augacc_y z_mysim*", full.names = T)
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
myresults_aipw_yz$Method = "acc-lambda-eta1" 

## combine both models
myresults_aipw = rbind(myresults_aipw_yz, myresults_aipw_z)
colnames(myresults_aipw) = c("beta1.x", "beta2.x", "beta0.x", "beta1.y", "beta2.y", "beta0.y", "Method")

head(myresults_aipw)

write.csv(myresults_aipw, "myresults2.csv", row.names = FALSE)


###############################################

myresults_aipw = rbind(read.csv("myresults1.csv"), read.csv("myresults2.csv"), read.csv("myresults3.csv"))


myresults_aipw = myresults_aipw %>%
  mutate(Method = factor(Method, levels = c("oracle", "naive", "cc", 
                                            "acc-lambda", "acc-lambda-eta1",  "acc-lambda-w", "acc-lambda-eta1-w"))) %>%
  mutate(Method.c = ifelse(Method %in% c("acc-lambda", "acc-lambda-eta1",  "acc-lambda-w", "acc-lambda-eta1-w"), "acc", 
                      ifelse(Method == "oracle", "oracle",
                      ifelse(Method == "cc", "cc", "naive")))) %>%
    mutate(Method.c = factor(Method.c, levels = c("oracle", "naive", "cc", "acc") ))
  
table(myresults_aipw$Method)
table(myresults_aipw$Method.c)

# beta0
p0 <- myresults_aipw %>%
  ggplot(aes(x=beta0.x, y=Method, fill=Method.c)) +
  geom_boxplot(outlier.color = "black", outlier.size = 0.2, position=position_dodge(.9)) +
  stat_summary(fun=mean,
               colour="white", geom="point", shape=18, size=2, position=position_dodge(.9)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() + theme(legend.position="bottom") +
  theme(axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.text.align = 0,
        axis.text.x=element_text(size=12)) +
  guides(fill=guide_legend(nrow=1)) +
  labs(y = "", x  = TeX("$\\widehat{\\beta}_0 = 1$"), fill = "Method") +
  scale_fill_manual(labels =c("Oracle", "Naive", "Complete Case", "ACC"), 
                    values=c("#000000",
                             "#E69F00",
                             #"#56B4E9",
                             #"#009E73",
                             "#F0E442",
                             "#0072B2",
                             "#D55E00"))  # "orange"
p0

ggsave(filename = paste0("1000-ACC-yz-beta0.png"), plot = p0,  width = 10, height = 6, dpi = 300, units = "in")


## radar plot
mysig = 95
zval  = 1-(1-mysig/100)/2
theta = c(0,0.2,0.2)
myresults_aipw$beta0.coverage <- ((myresults_aipw$beta0.x + qnorm(zval)*myresults_aipw$beta0.y) >= theta[1]) &
  ((myresults_aipw$beta0.x - qnorm(zval)*myresults_aipw$beta0.y) <= theta[1])
myresults_aipw$beta1.coverage <- ((myresults_aipw$beta1.x + qnorm(zval)*myresults_aipw$beta1.y) >= theta[2]) &
  ((myresults_aipw$beta1.x - qnorm(zval)*myresults_aipw$beta1.y) <= theta[2])
myresults_aipw$beta2.coverage <- ((myresults_aipw$beta2.x + qnorm(zval)*myresults_aipw$beta2.y) >= theta[3]) &
  ((myresults_aipw$beta2.x - qnorm(zval)*myresults_aipw$beta2.y) <= theta[3])
  

myresults_coverage = myresults_aipw %>%
  # calculate coverage
  mutate(beta0.coverage = ifelse(beta0.coverage, 1, 0),
         beta1.coverage = ifelse(beta1.coverage, 1, 0),
         beta2.coverage = ifelse(beta2.coverage, 1, 0)) %>%
  # select only variables of interest
  dplyr::select(Method, beta0.coverage, beta1.coverage, beta2.coverage) %>%
  # calculate means
  group_by(Method) %>%
  summarise_all(function(x) round(100*mean(x),2))

#transpose
myresults_coverage_t = t(myresults_coverage[-1,])[-1,] %>% as.numeric() %>% matrix(nrow = 3)
myresults_coverage_t = as.data.frame(myresults_coverage_t) 
colnames(myresults_coverage_t) = myresults_coverage$Method[-1]

library("ggradar")
myradar <- function(x=1, color = "#000000") {
  ggradar(
    myresults_coverage_t[x,],
    values.radar = c("", "", ""),
    grid.min = 0, grid.mid =95, grid.max = 150,
    font.radar = "Century Gothic",
    grid.label.size = 6,
    # Polygons
    group.line.width = 2,
    group.point.size = 5,
    group.colours = color,
    # Background and grid lines
    background.circle.colour = "white",
    gridline.min.colour = "darkgrey",
    gridline.mid.colour = "red"
  )
}


myradar(x=1, color = "#000000")
ggsave(filename = paste0("1000-ACC-yz-beta0-coverage.png"),  width = 12, height = 12, dpi = 300, units = "in")