# Exercise 9.27 and 10.29 using MRW1992

rm(list=ls())
cat("\014") 

# Setting the working directory
setwd("C:/Users/DELL/Desktop/582 HW1")

# Libraries
library(readxl) # to read excel xlsx files

# Section 8.12 of Hansen 

data <- read.table("MRW1992.txt",header=TRUE)
N <- matrix(data$N,ncol=1)
lndY <- matrix(log(data$Y85)-log(data$Y60),ncol=1)
lnY60 <- matrix(log(data$Y60),ncol=1)
lnI <- matrix(log(data$invest/100),ncol=1)
lnG <- matrix(log(data$pop_growth/100+0.05),ncol=1)
lnS <- matrix(log(data$school/100),ncol=1)
xx <- as.matrix(cbind(lnY60,lnI,lnG,lnS,matrix(1,nrow(lndY),1)))
x <- xx[N==1,]
y <- lndY[N==1]
n <- nrow(x)
k <- ncol(x)

# Unrestricted regression
invx <-solve(t(x)%*%x)
b_ols <- solve((t(x)%*%x),(t(x)%*%y))

beta_ols = b_ols

e_ols <- rep((y-x%*%beta_ols),times=k)
xe_ols <- x*e_ols
V_ols <- (n/(n-k))*invx %*% (t(xe_ols) %*% xe_ols) %*% invx
se_ols <- sqrt(diag(V_ols))

print(beta_ols)
print(se_ols)

# Constrained regression
R <- c(0,1,1,1,0)
iR <- invx %*% R %*% solve(t(R) %*% invx %*% R) %*% t(R)
b_cls <- b_ols - iR %*% b_ols
e_cls <- rep((y-x %*% b_cls),times=k)
xe_cls <- x*e_cls
V_tilde <- (n/(n-k+1)) * invx %*% (t(xe_cls) %*% xe_cls) %*% invx
V_cls <- V_tilde - iR %*% V_tilde - V_tilde %*% t(iR) +iR %*% V_tilde %*% t(iR) 
se_cls <- sqrt(diag(V_cls))

print(b_cls)
print(se_cls)

# Efficient minimum distance
Vr <- V_ols%*%R%*%solve(t(R)%*%V_ols%*%R)%*%t(R)
b_emd <- b_ols - Vr%*%b_ols
e_emd <- rep((y-x%*%b_emd),times=k)
xe_emd <- x*e_emd
V2 <- (n/(n-k+1))*invx%*%(t(xe_emd)%*%xe_emd)%*%invx
V_emd <- V2 - V2%*%R%*%solve(t(R)%*%V2%*%R)%*%t(R)%*%V2
se_emd <- sqrt(diag(V_emd))



# Exercise 9.27 - Estimate the unrestricted model 
# and test the hypothesis that the three coefficients sum to zero.

R <- c(0, 1, 1, 1, 0)       # Restriction matrix
cvalue <- 0                 # R' beta = c 

Rbetacvalue <- t(R) %*% beta_ols - cvalue 

RVR <- t(R) %*% V_ols %*% R
Wald <- n*t(Rbetacvalue) %*% solve(RVR, Rbetacvalue)
p.Wald <- pchisq(Wald, df = 1, lower.tail = FALSE)

summary.Wald <- cbind(Wald, p.Wald)
colnames(summary.Wald) <- c("Wald statistic", "Wald p-value")

print(summary.Wald)



# Exercise 10.29 - Let theta = sum of beta2 + beta3 + beta4 

# (a) - Estimate the regression by unrestricted least squares and 
# report standard errors calculated by asymptotic, jackknife and the bootstrap

library(readxl)
library(lmtest)
library(sandwich)
library(tidyverse)

rm(list=ls())
cat("\014")
set.seed(42)


# (a) Estimate the regression by OLS and se from asymptotic, jackknife and bootstrap

# read in data

MRW1992.raw <- read.table("MRW1992.txt", header = TRUE)

MRW <- MRW1992.raw %>% 
  filter(N == 1) %>% 
  mutate(log.growth = log(Y85) - log(Y60),
         log.gdp60 = log(Y60),
         log.invest = log(invest/100),
         log.ngd = log(pop_growth/100 + 0.05),
         log.school = log(school/100))

n = nrow(MRW)

# OLS
OLS.reg <- lm(log.growth ~ log.gdp60 + log.invest + log.ngd + log.school,
              data = MRW) # OLS 

V.HC1 <- vcovHC(OLS.reg, 
                type = "HC1")  # HC vcov matrix
se.HC1 <- sqrt(diag(V.HC1))

coeftest(OLS.reg, vcov=V.HC1)


# jackknife

MRW.id <- MRW %>%
  mutate(id = row_number())

beta.jackknife = map_dfr(c(1:nrow(MRW.id)), function(x){
  reg.jk <- lm(log.growth ~ log.gdp60 + log.invest + log.ngd + log.school,
               data = MRW.id %>%
                 filter(id != x))
  
  return(reg.jk$coefficients)
}) %>%
  mutate(id = row_number()) 

beta.mean.jk <- beta.jackknife %>%
  select(-id) %>%
  summarise_all(list('mean' = mean)) %>%
  pivot_longer(everything()) %>%
  pull(value)

vcov.list.jk <- map(c(1:nrow(MRW.id)), function(x){
  beta.jk.dataframe <- beta.jackknife %>%
    filter(id == x) %>%
    select(-id) %>%
    pivot_longer(everything()) %>%
    pull(value)
  
  return((beta.jk.dataframe - beta.mean.jk) %*% t(beta.jk.dataframe - beta.mean.jk))
})

vcov.jackknife <- ((n-1)/n) * Reduce("+", vcov.list.jk)
se.jackknife <- sqrt(diag(vcov.jackknife))


# Bootstrap

B <- 1000

MRW.nested <- MRW.id %>%
  group_nest(id)

beta.boot <- map_dfr(c(1:B), function(x){
  reg.boot <- lm(log.growth ~ log.gdp60 + log.invest + log.ngd + log.school,
                 data = MRW.nested %>%
                   slice_sample(n = nrow(MRW.nested), 
                                replace = TRUE) %>% 
                   unnest(data)) 
  
  return(reg.boot$coefficients)
}) %>%
  mutate(id = row_number())

beta.mean.boot <- beta.boot %>%
  select(-id) %>%
  summarise_all(list('mean' = mean)) %>%
  pivot_longer(everything()) %>%
  pull(value)

vcov.list.boot <- map(c(1:nrow(MRW.id)), function(x){
  beta.boot.dataframe <- beta.boot %>%
    filter(id == x) %>%
    select(-id) %>%
    pivot_longer(everything()) %>%
    pull(value)
  
  return((beta.boot.dataframe - beta.mean.boot) %*% t(beta.boot.dataframe - beta.mean.boot))
})

vcov.boot <- (1/(n-1)) * Reduce("+", vcov.list.boot)
se.boot <- sqrt(diag(vcov.boot))

summary.10.29a <- rbind(OLS.reg$coefficients, se.HC1, se.jackknife, se.boot)
rownames(summary.10.29a) <- rbind("estimate", "se.asymptotic", "se.jackknife", "se.bootstrap")


# (b) - Estimate theta and report standard errors calculated by asymptotic, jackknife and the bootstrap

theta <- OLS.reg$coefficients[3] + OLS.reg$coefficients[4] + OLS.reg$coefficients[5] 
names(theta) <- "theta"

vcov.theta <- function(vcov.matrix){
  return(sqrt(sum(vcov.matrix[3:5, 3:5])))
}

se.asymptotic <- vcov.theta(V.HC1)
se.jackknife <- vcov.theta(vcov.jackknife)
se.bootstrap <- vcov.theta(vcov.boot)

summary10.29b <- c(theta, se.asymptotic, se.jackknife, se.bootstrap)
names(summary10.29b) <- c("theta", "se.theta.asymptotic", "se.theta.jackknife", "se.theta.bootstrap")

# (c) Report confidence intervals for theta using the percentile and BC methods.

alpha = 0.05
alphas <- c(alpha/2, 1- alpha/2)

# evaluate z_alpha for each alpha values
z.alphas <- qnorm(alphas)

theta.hat.boot <- beta.boot$log.invest + beta.boot$log.ngd + beta.boot$log.school

# evaluate z0.star
p.star <- mean(theta.hat.boot <= theta)
z0.star <- qnorm(p.star)

# calculate x(alpha) for each alpha values
x.alphas <- pnorm(z.alphas + 2*z0.star)

theta.boot.quantile <- theta.hat.boot

q.percentile <- quantile(theta.boot.quantile, probs = x.alphas)


# BCa method

theta.hat.loo <- beta.jackknife$log.invest + beta.jackknife$log.ngd + beta.jackknife$log.school
theta.bar.loo <- mean(theta.hat.loo)

diff.jack <- theta.bar.loo - theta.hat.loo
a.jack.num <- sum(diff.jack^3)
a.jack.den <- 6 * sum(diff.jack^2)^(3/2)
a.jack <- a.jack.num/a.jack.den

BCa.correction <- (z.alphas + z0.star) / (1-a.jack*(z.alphas + z0.star))
x.alphas.correction <- pnorm(z0.star + BCa.correction)

q.BCa <- quantile(theta.hat.boot, probs = x.alphas.correction)


