# ECON582 HW1 - Hansen 
# Exercise 9.26 and 10.28 using Nerlove1963

rm(list=ls())
cat("\014") 

setwd("C:/Users/DELL/Desktop/582 HW1")  # Setting the working directory
set.seed(42)                            # is the answer seed to everything

# Libraries
library(readxl)     # to read excel xlsx files

# Exercise 9.26

# Data
nerlove.raw <- read_xlsx("Nerlove1963.xlsx")
nerlove.log <- log(nerlove.raw)
colnames(nerlove.log) <- c("logTC", "logQ", "logPL", "logPK", "logPF")

# (a) - Unrestricted OLS

nerlove.OLS <- lm(logTC ~ logQ + logPL + logPK + logPF, data=nerlove.log, x=TRUE, y=TRUE)
summary(nerlove.OLS)

# variables from OLS

n.nerlove <- nrow(nerlove.log)
y.nerlove <- as.matrix(nerlove.OLS$y)
x.nerlove <- as.matrix(nerlove.OLS$x)
beta.OLS <- nerlove.OLS$coefficients
ehat.OLS <- nerlove.OLS$residuals

# OLS standard errors

xehat.OLS <- x.nerlove * ehat.OLS                         # X*ehat
omega.hat.OLS <- t(xehat.OLS) %*% xehat.OLS / n.nerlove   # Var(X*ehat) matrix for hetero
Qhatxx.OLS <- (t(x.nerlove) %*% x.nerlove) / n.nerlove    # (X'X)/n
Qhatxx.inv.OLS <- solve(Qhatxx.OLS)                       # default solve = giving inverse of a

vcov.homo.OLS <- n.nerlove * vcov(nerlove.OLS)
vcov.hetero.OLS <- Qhatxx.inv.OLS %*% omega.hat.OLS %*% Qhatxx.inv.OLS

se.homo.OLS <- sqrt(diag(vcov.homo.OLS)/n.nerlove)
se.hetero.OLS <- sqrt(diag(vcov.hetero.OLS)/n.nerlove)

summary.OLS <- rbind(beta.OLS, se.homo.OLS, se.hetero.OLS)
print(summary.OLS)

# Before we head to (b) and so on, create minimum distance estimator
# because CLS is essentially minimum distance estimator with a specific W matrix

#### Minimum distance estimation function
n = n.nerlove

md.est <- function(x, y, beta.hat, vcov.hat, W.hat, R, cvalue){
  RB <- t(R) %*% beta.hat 
  WR <- solve(W.hat, R) 
  RWR <- t(R) %*% solve(W.hat, R)
  beta.md <- beta.hat - WR %*% solve(RWR, RB - cvalue)
  
  K <- diag(nrow = length(beta.hat)) - WR %*% solve(RWR, t(R))
  vcov.md <- K %*% vcov.hat %*% t(K)
  se.md <- sqrt(diag(vcov.md)/n)
  
  diff.hat.md <- beta.hat - beta.md
  J.md <- n * t(diff.hat.md) %*% W.hat %*% diff.hat.md
  return(list("coefficients" = beta.md, "vcov" = vcov.md, "se" = se.md, "J.md" = J.md, "K" = K))
}

# (b) Meaning of the restriction = function log C is HD1 to PL, PK and PF


# (c) - Constrained Least Square 
# Note that Constrained Least Square estimation is simply to set
# W.hat = inverse of vcov.OLS = solve(vcov.OLS) = (vcov.OLS)^{-1}

R <- c(0, 0, 1, 1, 1)       # Restriction matrix
cvalue <- 1                 # R' beta = c 

# Homoskedastic CLS 

nerlove.homo.CLS <- md.est(x.nerlove, y.nerlove, beta.OLS, vcov.homo.OLS, W.hat = solve(vcov.homo.OLS), R, cvalue)
beta.homo.CLS <- as.matrix(t(nerlove.homo.CLS$coefficients))
colnames(beta.homo.CLS) <- c("(Intercept)", "logQ", "logPL", "logPK", "logPF")
rownames(beta.homo.CLS) <- "beta.homo.CLS"

se.homo.CLS <- as.matrix(t(nerlove.homo.CLS$se))
colnames(se.homo.CLS) <- c("(Intercept)", "logQ", "logPL", "logPK", "logPF")
rownames(se.homo.CLS) <- "se.homo.CLS"


# Heteroskedastic CLS 
# Note that W.hat is still now (vcov.hetero.OLS)^{-1} 

nerlove.hetero.CLS <- md.est(x.nerlove, y.nerlove, beta.OLS, vcov.hetero.OLS, W.hat = solve(vcov.hetero.OLS), R, cvalue)
beta.hetero.CLS <- as.matrix(t(nerlove.hetero.CLS$coefficients))
colnames(beta.hetero.CLS) <- c("(Intercept)", "logQ", "logPL", "logPK", "logPF")
rownames(beta.hetero.CLS) <- "beta.hetero.CLS"

se.hetero.CLS <- as.matrix(t(nerlove.hetero.CLS$se))
colnames(se.hetero.CLS) <- c("(Intercept)", "logQ", "logPL", "logPK", "logPF")
rownames(se.hetero.CLS) <- "se.hetero.CLS"

# CLS summary :
summary.CLS <- rbind(beta.homo.CLS, se.homo.CLS, beta.hetero.CLS, se.hetero.CLS)
print(summary.CLS)


# (d) Efficient Minimum Distance 

# Homoskedastic emd is the same as homoskedastic CLS 

nerlove.homo.emd <- md.est(x.nerlove, y.nerlove, beta.OLS, vcov.homo.OLS, W.hat = solve(vcov.homo.OLS), R, cvalue)
beta.homo.emd <- as.matrix(t(nerlove.homo.emd$coefficients))
colnames(beta.homo.emd) <- c("(Intercept)", "logQ", "logPL", "logPK", "logPF") 
rownames(beta.homo.emd) <- "beta.homo.emd" 

se.homo.emd <- as.matrix(t(nerlove.homo.emd$se))
colnames(se.homo.emd) <- c("(Intercept)", "logQ", "logPL", "logPK", "logPF") 
rownames(se.homo.emd) <- "se.homo.emd" 

# Heteroskedastic emd uses different W matrix = (vcov.hetero.OLS)^{-1}

nerlove.hetero.emd <- md.est(x.nerlove, y.nerlove, beta.OLS, vcov.hetero.OLS, W.hat = solve(vcov.hetero.OLS), R, cvalue)
beta.hetero.emd <- as.matrix(t(nerlove.hetero.emd$coefficients))
colnames(beta.hetero.emd) <- c("(Intercept)", "logQ", "logPL", "logPK", "logPF")
rownames(beta.hetero.emd) <- "beta.hetero.emd"

se.hetero.emd <- as.matrix(t(nerlove.hetero.emd$se))
colnames(se.hetero.emd) <- c("(Intercept)", "logQ", "logPL", "logPK", "logPF")
rownames(se.hetero.emd) <- "se.hetero.emd"

# emd summary :
summary.emd <- rbind(beta.homo.emd, se.homo.emd, beta.hetero.emd, se.hetero.emd)
print(summary.emd)


# (e) Testing the restriction using Wald statistic

Rbetacvalue <- t(R) %*% beta.OLS - cvalue 
RVR.homo <- t(R) %*% vcov.homo.OLS %*% R
Wald.homo <- n * t(Rbetacvalue) %*% solve(RVR.homo, Rbetacvalue)
p.Wald.homo <- pchisq(Wald.homo, df = 1, lower.tail = FALSE)

RVR.hetero <- t(R) %*% vcov.hetero.OLS %*% R
Wald.hetero <- n * t(Rbetacvalue) %*% solve(RVR.hetero, Rbetacvalue)
p.Wald.hetero <- pchisq(Wald.hetero, df = 1, lower.tail = FALSE)

summary.Wald <- rbind(cbind(Wald.homo, p.Wald.homo), cbind(Wald.hetero, p.Wald.hetero))
colnames(summary.Wald) <- c("Wald statistic", "p-value")
rownames(summary.Wald) <- rbind("under Homo", "under Hetero")

print(summary.Wald) 

# (f) Testing the restriction using Minimum Distance statistic 

J.CLS.homo <- nerlove.homo.CLS$J.md
p.J.CLS.homo <- pchisq(J.CLS.homo, df=1, lower.tail = FALSE)

J.CLS.hetero <- nerlove.hetero.CLS$J.md
p.J.CLS.hetero <- pchisq(J.CLS.hetero, df=1, lower.tail = FALSE)

# Note that J.CLS.homo = hetero because W.hat matrix does not change in CLS 

J.emd.homo <- nerlove.homo.emd$J.md
p.J.emd.homo <- pchisq(J.emd.homo, df=1, lower.tail = FALSE)

J.emd.hetero <- nerlove.hetero.emd$J.md
p.J.emd.hetero <- pchisq(J.emd.hetero, df=1, lower.tail = FALSE)


summary.md <- rbind(cbind(J.CLS.homo, p.J.CLS.homo, J.CLS.hetero, p.J.CLS.hetero), 
                   cbind(J.emd.homo, p.J.emd.homo, J.emd.hetero, p.J.emd.hetero))
colnames(summary.md) <- c("md Statistic under homo", "p-value under homo", "md Statistic under hetero", "p-value under hetero")
rownames(summary.md) <- rbind("using CLS", "using emd")

print(summary.md)



# Exercise 10.28

# (a) Estimate the regression by OLS and se from asymptotic, jackknife and bootstrap

# Theta = beta3 + beta4 + beta5

theta.est <- function(data){
  data.OLS <- lm(logTC ~ logQ + logPL + logPK + logPF, data = data)
  
  beta.hat <- coef(data.OLS)
  sig2.hat <- mean(resid(data.OLS)^2)
  theta.hat <-
    beta.hat[3] + beta.hat[4] + beta.hat[5]
  
  return(theta.hat)
}

# Jackknife

data_id <- nerlove.log %>%
  mutate(id = row_number()) # make id column

theta.hat.loo <- rep(0, n)

for(i in 1:n){
  df.loo.i <- data_id[-i, ]
  theta.hat.loo[i] <- theta.est(data = df.loo.i)
}

theta.bar.loo <- mean(theta.hat.loo) 
var.jackknife <- (n-1) * mean((theta.hat.loo - theta.bar.loo)^2)
se.jackknife <- sqrt(var.jackknife)



# (b) Estimate theta = beta3 + beta4 + beta5 and se from asymptotic, jackknife and bootstrap


# (c) Report confidence intervals for theta using the percentile and BCa methods









