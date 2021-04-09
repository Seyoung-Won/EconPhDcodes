# ECON582 HW1 - Hansen 
# Exercise 4.26 and 10.31 using DDK2011

# Setting the working directory

setwd("C:/Users/DELL/Desktop/582 HW1")
set.seed(42)


# Libraries
library(readxl) # to read excel xlsx files
library(lmtest)
library(sandwich)
library(tidyverse)
library(data.table)

# Clear the environment and console

rm(list=ls())
cat("\014") 

# Hansen 4.26

# Data
DDK.raw <- read_excel("DDK2011.xlsx", na=".")

DDK4.26 <- DDK.raw %>%
  mutate(totalscore_z = ((totalscore - mean(totalscore)) / sd(totalscore))) %>%
  select(totalscore_z, tracking, agetest, girl, etpteacher, schoolid, percentile) %>%
  mutate_all(as.numeric) 

DDK4.26 <- DDK.raw %>%
  na.omit(data.frame()) %>%
  mutate(totalscore_z = ((totalscore - mean(totalscore)) / sd(totalscore))) %>%
  select(totalscore_z, tracking, agetest, girl, etpteacher, schoolid, percentile) %>%
  mutate_all(as.numeric) 


basic.ols <- lm(totalscore_z ~ tracking + agetest + girl + etpteacher + percentile, 
                data = DDK4.26)

vcov.HC1 <- vcovHC(basic.ols, 
                type = "HC1")
se.HC1 <- sqrt(diag(vcov.HC1))

vcov.cluster <- vcovCL(basic.ols,
                       cluster = ~ schoolid)
se.cluster <- sqrt(diag(vcov.cluster))

summary.4.26 <- cbind(basic.ols$coefficients, 
                      se.HC1, 
                      se.cluster,
                      se.HC1 - se.cluster)
  colnames(summary.4.26) <- c("estimate", "se.HC1", "se.cluster", "se.difference")



# Hansen 10.31 

# (a) bootstrap
  
DDK10.31 <- DDK4.26
DDK.nested <- DDK4.26 %>%
  group_nest(schoolid)

B <- 1000

beta.boot <- map_dfr(c(1:B), function(data){
  ols.reg <- lm(totalscore_z ~ tracking + agetest + girl + etpteacher + percentile,
                data = DDK.nested %>% 
                  slice_sample(n = nrow(DDK.nested), replace = TRUE) %>%
                  unnest(data))
  
  return(ols.reg$coefficients)
})

mean.boot <- beta.boot %>%
  summarise_all(list('mean' = mean)) %>%
  pivot_longer(everything()) %>%
  pull(value)

vcov.list <- map(c(1:nrow(beta.boot)), function(x){
  df <- beta.boot %>%
    filter(row_number() == x) %>%
    pivot_longer(everything()) %>%
    pull(value)
  
  return((df - mean.boot) %*% t(df - mean.boot))
})

var.boot <- (1/(B-1)) * Reduce("+", vcov.list)
se.boot <- sqrt(diag(var.boot))

summary.10.31 <- cbind(basic.ols$coefficients, 
                        se.HC1, 
                        se.cluster,
                        se.boot)
colnames(summary.10.31) <- c("estimate", "se.HC1", "se.cluster", "se.bootstrap")

summary.10.31

# BCa for all coefficients

n = nrow(DDK10.31)

beta.jackknife <- map_dfr(c(1:n), function(x){
  reg.jk <- lm(totalscore_z ~ tracking + agetest + girl + etpteacher + percentile,
                  data = DDK10.31 %>%
                  filter(row_number() != x))
  
  return(reg.jk$coefficients)
}) %>%
  mutate(id = row_number())

# CI 

alpha = 0.05
alphas <- c(alpha/2, 1- alpha/2)

# evaluate z_alpha for each alpha values
z.alphas <- qnorm(alphas)

# evaluate z0.star
beta.hat <- basic.ols$coefficients

beta.boot.CI <- beta.boot %>%
  as.matrix

p.star.intercept <- mean(beta.boot[1] <= beta.hat[1])
p.star.tracking <- mean(beta.boot[2] <= beta.hat[2])
p.star.agetest <- mean(beta.boot[3] <= beta.hat[3])
p.star.girl <- mean(beta.boot[4] <= beta.hat[4])
p.star.etpteacher <- mean(beta.boot[5] <= beta.hat[5])
p.star.percentile <- mean(beta.boot[6] <= beta.hat[6])

p.stars <- c(p.star.intercept, p.star.tracking, p.star.agetest, 
             p.star.girl, p.star.etpteacher, p.star.percentile)

z0.stars <- qnorm(p.stars)

# calculate x(alpha) for each alpha values

CI <- map_dfc(c(1:6), function(x){
  
  x.alphas <- pnorm(z.alphas + 2*z0.stars[x])
  
  beta.jk <- beta.jackknife %>%
    as.matrix
  
  diff.jack <- beta.jk[,x] - mean(beta.jk[,x])
  a.jack.num <- sum(diff.jack^3)
  a.jack.den <- 6 * sum(diff.jack^2)^(3/2)
  a.jack <- a.jack.num/a.jack.den
  
  BCa.correction <- (z.alphas + z0.stars[x]) / (1-a.jack*(z.alphas + z0.stars[x]))
  x.alphas.correction <- pnorm(z0.stars[x] + BCa.correction)
  
  q.BCa <- quantile(beta.boot.CI[,x], probs = x.alphas.correction) %>%
    as.data.frame %>%
    rownames_to_column(var = "quantile")
  
  return(q.BCa)
})
colnames(CI) <- c("Intercept Percentile", "Intercept bounds",
                  "tracking Percentile", "tracking bounds",
                  "agetest Percentile", "agetest bounds",
                  "girl Percentile", "girl bounds",
                  "etpteacher Percentile", "etpteacher bounds",
                  "percentile Percentile", "percentile bounds")

print(CI)

