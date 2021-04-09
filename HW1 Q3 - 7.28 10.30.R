# Exercise 7.28 and 10.30 using cps09mar

# Setting the working directory and random seed
setwd("C:/Users/DELL/Desktop/582 HW1")
set.seed(42) 

# Libraries
library(readxl)     # to read excel xlsx files
library(sandwich)   # make nice HC covariance matrices
library(lmtest)     # for robust coefficient testing
library(tidyverse)  # important for data manipulation, includes dplyr

# Exercise 7.28

rm(list=ls())
cat("\014") 

# Data
cps09mar.raw <- read_xlsx("cps09mar.xlsx")

cps09mar <- cps09mar.raw %>%
  mutate(log.wage = log(earnings/(hours*week))) %>%
  mutate(experience = age-education-6) %>%
  mutate(experience2.100 = experience^2/100) %>%
  filter(race==1, female==0, hisp==1)


# (a) Report the coefficient estimates and robust standard errors.

OLS.reg <- lm(log.wage ~ education + experience + experience2.100,
              data = cps09mar) # OLS 

V.HC1 <- vcovHC(OLS.reg, 
                type = "HC1")  # HC vcov matrix

coeftest(OLS.reg, vcov=V.HC1)


# (b) theta = beta1 (beta2 + beta3/5)

b1 <- coef(OLS.reg)[2]
b2 <- coef(OLS.reg)[3]
b3 <- coef(OLS.reg)[4]

theta.hat = b1 / (b2  + b3/5)
names(theta.hat) <- "theta.hat"

# (c) se of theta as a function of V.HC1 and se.theta

# Restriction vector = p.derivatives of theta w.r.t b0 b1 b2 b3
R <- c(0, 1/(b2 + b3/5), -b1/(b2+b3/5)^2, (-1/5)*b1/(b2+b3/5)^2)

# V.theta = scalar
vcov.theta <- t(R) %*% V.HC1 %*% R %>%
  as.numeric

se.theta <- sqrt(vcov.theta)

# (d) construct 90% CI

CI90 <- c(theta.hat - qnorm(0.95)*se.theta, theta.hat + qnorm(0.95)*se.theta)
names(CI90) <- c("5%", "95%")


# (e) regression at edu=12 exp=20 and 95% CI
# We already have the beta estimates and V.HC1. So just plug-in

edu12 <- 12
exp20 <- 20

plugin.vector <- rbind(1, edu12, exp20, exp20^2/100)

# y = X'beta, here beta = column vector
plugin.value <- t(plugin.vector) %*% (OLS.reg$coefficients) %>%
  as.numeric

# getting the scalar se
vcov.plugin <- t(plugin.vector) %*% V.HC1 %*% (plugin.vector)
se.plugin <- sqrt(vcov.plugin) %>%
  as.numeric

CI95 <- c(plugin.value - qnorm(0.975)*se.plugin, plugin.value + qnorm(0.975)*se.plugin)
names(CI95) <- c("2.5%", "97.5%")


# (f) Consider out-of-sample individuals with edu=16 exp=5 
#     Construct 80% forecast interval for their log wage and wage

edu16 <- 16
exp5 <- 5

forecast.vector <- rbind(1, edu16, exp5, exp5^2/100)

forecast.value <- t(forecast.vector) %*% (OLS.reg$coefficients) %>%
  as.numeric

# Forecast intervals reguires sigma.hat2 (Hansen 7.15)

sigma2 <- mean((OLS.reg$residuals)^2) # (E(residual))^2 = 0

vcov.forecast <- sigma2 + t(forecast.vector) %*% V.HC1 %*% (forecast.vector)
se.forecast <- sqrt(vcov.forecast) %>%
  as.numeric

CI80 <- c(forecast.value - qnorm(0.9)*se.forecast, forecast.value + qnorm(0.9)*se.forecast)
names(CI80) <- c("10%", "90%")

CI80exp <- c(exp(forecast.value - qnorm(0.9)*se.forecast), exp(forecast.value + qnorm(0.9)*se.forecast))
names(CI80exp) <- c("10%", "90%")


# Exercise 10.30

rm(list=ls())
cat("\014") 

# (a) - Estimate theta and report standard errors calculated 
# by asymptotic, jackknife and the bootstrap.

cps09mar.raw <- read_xlsx("cps09mar.xlsx")

# For MidWest region
cps09mar.MW <- cps09mar.raw %>%
  mutate(log.wage = log(earnings/(hours*week))) %>%
  mutate(experience = age-education-6) %>%
  mutate(experience2.100 = experience^2/100) %>%
  filter(race==1, female==0, hisp==1) %>%
  filter(marital==7, region==2) 

OLS.reg.MW <- lm(log.wage ~ education + experience + experience2.100,
              data = cps09mar.MW) # OLS 

V.HC1.MW <- vcovHC(OLS.reg.MW, 
                type = "HC1")  # HC vcov matrix

coeftest(OLS.reg.MW, vcov=V.HC1.MW)

b1.MW <- coef(OLS.reg.MW)[2]
b2.MW <- coef(OLS.reg.MW)[3]
b3.MW <- coef(OLS.reg.MW)[4]

theta.hat.MW = b1.MW / (b2.MW  + b3.MW/5)
names(theta.hat.MW) <- "theta.hat.Midwest"

# Asymptotic 

  # Restriction vector = p.derivatives of theta w.r.t b0 b1 b2 b3
  R.MW <- c(0, 1/(b2.MW + b3.MW/5), -b1.MW/(b2.MW+b3.MW/5)^2, (-1/5)*b1.MW/(b2.MW+b3.MW/5)^2)
  
  # V.theta = scalar
  vcov.theta.asym <- t(R.MW) %*% V.HC1.MW %*% R.MW %>%
    as.numeric

  se.theta.asym <- sqrt(vcov.theta.asym)
  names(se.theta.asym) <- "se.theta.asymptotic"

# Jackknife

n = nrow(cps09mar.MW)
  
cps09mar.MW.id <- cps09mar.MW %>%
  mutate(id = row_number())
  
theta.jk = map_dfr(c(1:n), function(x){
  reg.jk <- lm(log.wage ~ education + experience + experience2.100,
               data = cps09mar.MW.id %>%
                 filter(id != x))
  
  b1.jk <- reg.jk$coefficient[2]
  b2.jk <- reg.jk$coefficient[3]
  b3.jk <- reg.jk$coefficient[4]
  
  theta.hat.jk <- b1.jk / (b2.jk + b3.jk/5)
  names(theta.hat.jk) <- "theta.hat.jk"
  
  return(theta.hat.jk)
}) %>%
  mutate(id = row_number())

theta.mean.jk <- theta.jk %>%
  select(-id) %>%
  summarise_all(list('mean' = mean)) %>%
  pivot_longer(everything()) %>%
  pull(value)
  
vcov.list.jk <- map(c(1:n), function(x){
  theta.jk.dataframe <- theta.jk %>%
    filter(id == x) %>%
    select(-id) %>%
    pivot_longer(everything()) %>%
    pull(value)
  
  return((theta.jk.dataframe - theta.mean.jk) %*% t(theta.jk.dataframe - theta.mean.jk))
})
  
vcov.theta.jk <- ((n-1)/n) * Reduce("+", vcov.list.jk)
se.theta.jackknife <- sqrt(diag(vcov.theta.jk))
names(se.theta.jackknife) <- "se.theta.jackknife"
  

# Bootstrap

B <- 1000
  
cps09mar.MW.id.nested <- cps09mar.MW.id %>%
  group_nest(id)

theta.boot <- map_dfr(c(1:B), function(x){
  reg.boot <- lm(log.wage ~ education + experience + experience2.100,
               data = cps09mar.MW.id.nested %>%
                 slice_sample(n = nrow(cps09mar.MW.id.nested),
                              replace = TRUE) %>%
                 unnest(data))
  
  b1.boot <- reg.boot$coefficient[2]
  b2.boot <- reg.boot$coefficient[3]
  b3.boot <- reg.boot$coefficient[4]
  
  theta.hat.boot <- b1.boot / (b2.boot + b3.boot/5)
  names(theta.hat.boot) <- "theta.hat.boot"
  
  return(theta.hat.boot)
}) %>%
  mutate(id = row_number())

theta.mean.boot <- theta.boot %>%
  select(-id) %>%
  summarise_all(list('mean' = mean)) %>%
  pivot_longer(everything()) %>%
  pull(value)

vcov.list.boot <- map(c(1:n), function(x){
  theta.boot.dataframe <- theta.boot %>%
    filter(id == x) %>%
    select(-id) %>%
    pivot_longer(everything()) %>%
    pull(value)
  
  return((theta.boot.dataframe - theta.mean.boot) %*% t(theta.boot.dataframe - theta.mean.boot))
})

vcov.boot <- (1/(n-1)) * Reduce("+", vcov.list.boot)
se.theta.boot <- sqrt(diag(vcov.boot))
names(se.theta.boot) <- "se.theta.bootstrap"

summary.10.30a <- c(theta.hat.MW, se.theta.asym, se.theta.jackknife, se.theta.boot)



# (b) - Explain the discrepancy between the standard errors.

# Sample size only n = 99 which is very small, compared to n = 4230 in Ex. 7.28
# This sample size may not be large enough to have LLN to kick in at asymptote.
# Bootstrap has a greater se because bootstrap replaces the samples so it has a higher prob. of outliers


# (c) - Report confidence intervals for theta using the BC percentile method

alpha = 0.05
alphas <- c(alpha/2, 1- alpha/2)

# evaluate z_alpha for each alpha values
z.alphas <- qnorm(alphas)

# evaluate z0.star
p.star <- mean(theta.boot <= theta.hat.MW)
z0.star <- qnorm(p.star)

# calculate x(alpha) for each alpha values
x.alphas <- pnorm(z.alphas + 2*z0.star)

theta.boot.quantile <- theta.boot %>%
  select(-id) %>%
  as.matrix

quantile(theta.boot.quantile, probs = x.alphas)

