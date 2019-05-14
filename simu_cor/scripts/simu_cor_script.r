# simulation analysis of trade-off / synergies, aka correlation / covariation 

## set wd
setwd("~/Desktop/Blog/GitFolder/simu_cor/")

## load libraries
library(MASS) # for simulation
library(brms) # for modelling
library(plyr)
library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


## helper function

# a function to generate covariance matrices from correlation matrix and stdev
cor2cov <- function(cor, sds){
  # code from: https://stackoverflow.com/questions/18740796/generate-covariance-matrix-from-correlation-matrix
  cov <- R * tcrossprod(sds)
  return(cov)
}

# a function to grab the an estimated variance-covariance matrix from a brms object
cov_brms <- function(model, fun = "median"){
  post <- posterior_samples(model)
  # grab estimate of sd
  sd_var <- apply(post[,grep("sigma",names(post))], 2, match.fun(fun))
  # grab estimate of correlation
  r_var <- apply(post[,grep("rescor",names(post))], 2, match.fun(fun))
  # put this in a matrix format
  R <- data.frame(X = names(r_var),Cor = r_var)
  R <- separate(R, X, c("drop","X1","X2"))
 
  
  R2 <- dcast(R,X1~X2,value.var="Cor")
  rownames(R2) <- R2[,1]
  R2 <- R2[,-1]
  R2 <- cbind(NA, R2)
  #colnames(R2)[1] <- "Cstock"
  R2 <- rbind(R2,NA)
  #rownames(R2)[nrow(R2)] <- "Predation"
  R2 <- as.matrix(R2)
  R2[lower.tri(R2)] <- t(R2)[lower.tri(R2)]
  diag(R2) <- sd_var
  
  return(R2)
}

## first simulation, two response variables, one predictor with positive effect on both and positive covariation
# the (residual) correlation matrix
R <- matrix(c(1, 0.6, 0.6, 1),ncol=2)
# the variable sd
stdev <- c(2, 1.5)
# get the residual covariance matrix
D <- cor2cov(R, stdev)
# the predictor
X <- runif(100, -2, 2)
# the linear predictor with different intercept for the two variables but same slope
mu <- cbind(1 + 2 * X, 5 + 2 * X)
# simulate the data packed in a function
sim_exp1 <- function(){
  Y <- adply(mu, 1, function(mu) rmvnorm(n = 1, mean = mu, sigma = D))
  # the observed correlation matrix
  bias_obs <- data.frame(Bias = cor(Y[,-1])[2,1] - R[2,1], Type = "observed")
  
  # format for brms
  brm_dat <- data.frame(X = X, Y1 = Y[,2], Y2 = Y[,3])
  # fit the model
  to_stan <- make_standata(mvbind(Y1, Y2) ~ X, data = brm_dat)
  m_mult <- stan("models/mult_model.stan",data = to_stan)
  # get the fitted residual correlation
  bias_m <- data.frame(Bias = median(extract(m_mult)$rescor - R[2,1]), Type = "model")
  bias_all <- rbind(bias_obs, bias_m)
  return(bias_all)
}

# replicate this 100 times
rep_exp1 <- replicate(100,sim_exp1(), simplify = FALSE)
rep_exp1 <- rbind.fill(rep_exp1)
rep_exp1$real_R <- 0.6
rep_exp1$fit_R <- with(rep_exp1, Bias + real_R)
rep_exp1$Scenario <- "Scenario 1"

ggplot(rep_exp1,aes(x=Type,y=fit_R)) +
  geom_jitter() +
  geom_hline(yintercept = rep_exp1$real_R, linetype = "dashed")

## second simulation, two response variables, one predictor, with positive effect on both and null residual correlation
# the (residual) correlation matrix
R <- matrix(c(1, 0.01, 0.01, 1),ncol=2)
# the variable sd
stdev <- c(2, 1.5)
# get the residual covariance matrix
D <- cor2cov(R, stdev)
# the predictor
X <- runif(100, -2, 2)
# the linear predictor with different intercept for the two variables but same slope
mu <- cbind(1 + 2 * X, 5 + 2 * X)
# simulate the data packed in a function
sim_exp2 <- function(){
  Y <- adply(mu, 1, function(mu) rmvnorm(n = 1, mean = mu, sigma = D))
  # the observed correlation matrix
  bias_obs <- data.frame(Bias = cor(Y[,-1])[2,1] - R[2,1], Type = "observed")
  
  # format for brms
  brm_dat <- data.frame(X = X, Y1 = Y[,2], Y2 = Y[,3])
  # fit the model
  to_stan <- make_standata(mvbind(Y1, Y2) ~ X, data = brm_dat)
  m_mult <- stan("models/mult_model.stan",data = to_stan)
  # get the fitted residual correlation
  bias_m <- data.frame(Bias = median(extract(m_mult)$rescor - R[2,1]), Type = "model")
  bias_all <- rbind(bias_obs, bias_m)
  return(bias_all)
}

# replicate this 100 times
rep_exp2 <- replicate(100,sim_exp2(),simplify = FALSE)
rep_exp2 <- rbind.fill(rep_exp2)
rep_exp2$real_R <- 0.01
rep_exp2$fit_R <- with(rep_exp2, Bias + real_R)
rep_exp2$Scenario <- "Scenario 2"

ggplot(rep_exp2,aes(x=Type,y=fit_R)) +
  geom_jitter() +
  geom_hline(yintercept = rep_exp2$real_R, linetype = "dashed")

## third simulation, two response variables, one predictor, with conteracting effect and positive residual correlation
# the (residual) correlation matrix
R <- matrix(c(1, 0.6, 0.6, 1),ncol=2)
# the variable sd
stdev <- c(2, 1.5)
# get the residual covariance matrix
D <- cor2cov(R, stdev)
# the predictor
X <- runif(100, -2, 2)
# the linear predictor with different intercept for the two variables but same slope
mu <- cbind(1 + 0.5 * X, 5 - 0.5 * X)
# simulate the data packed in a function
sim_exp3 <- function(){
  Y <- adply(mu, 1, function(mu) rmvnorm(n = 1, mean = mu, sigma = D))
  # the observed correlation matrix
  bias_obs <- data.frame(Bias = cor(Y[,-1])[2,1] - R[2,1], Type = "observed")
  
  # format for brms
  brm_dat <- data.frame(X = X, Y1 = Y[,2], Y2 = Y[,3])
  # fit the model
  to_stan <- make_standata(mvbind(Y1, Y2) ~ X, data = brm_dat)
  m_mult <- stan("models/mult_model.stan",data = to_stan)
  # get the fitted residual correlation
  bias_m <- data.frame(Bias = median(extract(m_mult)$rescor - R[2,1]), Type = "model")
  bias_all <- rbind(bias_obs, bias_m)
  return(bias_all)
}

# replicate this 100 times
rep_exp3 <- replicate(100,sim_exp3(),simplify=FALSE)
rep_exp3 <- rbind.fill(rep_exp3)
rep_exp3$real_R <- 0.6
rep_exp3$fit_R <- with(rep_exp3, Bias + real_R)
rep_exp3$Scenario <- "Scenario 3"

ggplot(rep_exp3,aes(x=Type,y=fit_R)) +
  geom_jitter() +
  geom_hline(yintercept = rep_exp3$real_R, linetype = "dashed")

## fourth simulation, two response variables, one predictor, with conteracting effect and no residual correlation
# the (residual) correlation matrix
R <- matrix(c(1, -0.01, -0.01, 1),ncol=2)
# the variable sd
stdev <- c(2, 1.5)
# get the residual covariance matrix
D <- cor2cov(R, stdev)
# the predictor
X <- runif(100, -2, 2)
# the linear predictor with different intercept for the two variables but same slope
mu <- cbind(1 + 2 * X, 5 - 2 * X)
# simulate the data packed in a function
sim_exp4 <- function(){
  Y <- adply(mu, 1, function(mu) rmvnorm(n = 1, mean = mu, sigma = D))
  # the observed correlation matrix
  bias_obs <- data.frame(Bias = cor(Y[,-1])[2,1] - R[2,1], Type = "observed")
  
  # format for brms
  brm_dat <- data.frame(X = X, Y1 = Y[,2], Y2 = Y[,3])
  # fit the model
  to_stan <- make_standata(mvbind(Y1, Y2) ~ X, data = brm_dat)
  m_mult <- stan("models/mult_model.stan",data = to_stan)
  # get the fitted residual correlation
  bias_m <- data.frame(Bias = median(extract(m_mult)$rescor - R[2,1]), Type = "model")
  bias_all <- rbind(bias_obs, bias_m)
  return(bias_all)
}

# replicate this 100 times
rep_exp4 <- replicate(100,sim_exp4(),simplify = FALSE)
rep_exp4 <- rbind.fill(rep_exp4)
rep_exp4$real_R <- -0.01
rep_exp4$fit_R <- with(rep_exp4, Bias + real_R)
rep_exp4$Scenario <- "Scenario 4"

ggplot(rep_exp4,aes(x=Type,y=fit_R)) +
  geom_jitter() +
  geom_hline(yintercept = rep_exp4$real_R, linetype = "dashed")

# plot of these first 4 simulations
rep_a <- rbind(rep_exp1, rep_exp2, rep_exp3, rep_exp4)
rep_a %>%
  group_by(Scenario) %>%
  summarise(real_R = mean(real_R)) -> rep_dd

gg_1 <- ggplot(rep_a,aes(x=Type,y=fit_R)) +
  geom_jitter() +
  geom_hline(aes(yintercept = real_R), linetype = "dashed") +
  facet_wrap(~Scenario)
  
ggsave("figures/simulation1_4_res.png",gg_1)


## fifth simulation, two response variables, two predictors, with same effect on the responses and no residual correlation, two predictors in the model
# the (residual) correlation matrix
R <- matrix(c(1, -0.01, -0.01, 1),ncol=2)
# the variable sd
stdev <- c(2, 1.5)
# get the residual covariance matrix
D <- cor2cov(R, stdev)
# the predictor
X <- runif(100, -2, 2)
Z <- runif(100, -2, 2)
# the linear predictor with different intercept for the two variables but same slope
mu <- cbind(1 + 2 * X + Z, 5 + 2 * X + Z)
# simulate the data packed in a function
sim_exp5 <- function(){
  Y <- adply(mu, 1, function(mu) rmvnorm(n = 1, mean = mu, sigma = D))
  # the observed correlation matrix
  bias_obs <- data.frame(Bias = cor(Y[,-1])[2,1] - R[2,1], Type = "observed")
  
  # format for brms
  brm_dat <- data.frame(X = X, Z = Z, Y1 = Y[,2], Y2 = Y[,3])
  # fit the model
  to_stan <- make_standata(mvbind(Y1, Y2) ~ X + Z, data = brm_dat)
  m_mult <- stan("models/mult_model.stan",data = to_stan)
  # get the fitted residual correlation
  bias_m <- data.frame(Bias = median(extract(m_mult)$rescor - R[2,1]), Type = "model")
  bias_all <- rbind(bias_obs, bias_m)
  return(bias_all)
}

# replicate this 100 times
rep_exp5 <- replicate(100,sim_exp5(),simplify = FALSE)
rep_exp5 <- rbind.fill(rep_exp5)
rep_exp5$real_R <- -0.01
rep_exp5$fit_R <- with(rep_exp5, Bias + real_R)


ggplot(rep_exp5,aes(x=Type,y=fit_R)) +
  geom_jitter() +
  geom_hline(yintercept = rep_exp5$real_R, linetype = "dashed")


## sixth simulation, two response variables, two predictors, with differing effect on the responses and no residual correlation, two predictors in the model
# the (residual) correlation matrix
R <- matrix(c(1, -0.01, -0.01, 1),ncol=2)
# the variable sd
stdev <- c(2, 1.5)
# get the residual covariance matrix
D <- cor2cov(R, stdev)
# the predictor
X <- runif(100, -2, 2)
Z <- runif(100, -2, 2)
# the linear predictor with different intercept for the two variables but same slope
mu <- cbind(1 + 2 * X + Z, 5 - 2 * X - Z)
# simulate the data packed in a function
sim_exp6 <- function(){
  Y <- adply(mu, 1, function(mu) rmvnorm(n = 1, mean = mu, sigma = D))
  # the observed correlation matrix
  bias_obs <- data.frame(Bias = cor(Y[,-1])[2,1] - R[2,1], Type = "observed")
  
  # format for brms
  brm_dat <- data.frame(X = X, Z = Z, Y1 = Y[,2], Y2 = Y[,3])
  # fit the model
  to_stan <- make_standata(mvbind(Y1, Y2) ~ X + Z, data = brm_dat)
  m_mult <- stan("models/mult_model.stan",data = to_stan)
  # get the fitted residual correlation
  bias_m <- data.frame(Bias = median(extract(m_mult)$rescor - R[2,1]), Type = "model")
  bias_all <- rbind(bias_obs, bias_m)
  return(bias_all)
}

# replicate this 100 times
rep_exp6 <- replicate(100,sim_exp6(),simplify = FALSE)
rep_exp6 <- rbind.fill(rep_exp6)
rep_exp6$real_R <- -0.01
rep_exp6$fit_R <- with(rep_exp6, Bias + real_R)

ggplot(rep_exp6,aes(x=Type,y=fit_R)) +
  geom_jitter() +
  geom_hline(yintercept = rep_exp6$real_R, linetype = "dashed")

## seventh simulation, two response variables, two predictors, with differing effect on the responses and no residual correlation, one predictors in the model
# the (residual) correlation matrix
R <- matrix(c(1, -0.01, -0.01, 1),ncol=2)
# the variable sd
stdev <- c(2, 1.5)
# get the residual covariance matrix
D <- cor2cov(R, stdev)
# the predictor
X <- runif(100, -2, 2)
Z <- runif(100, -2, 2)
# the linear predictor with different intercept for the two variables but same slope
mu <- cbind(1 + 2 * X + Z, 5 - 2 * X - Z)
# simulate the data packed in a function
sim_exp7 <- function(){
  Y <- adply(mu, 1, function(mu) rmvnorm(n = 1, mean = mu, sigma = D))
  # the observed correlation matrix
  bias_obs <- data.frame(Bias = cor(Y[,-1])[2,1] - R[2,1], Type = "observed")
  
  # format for brms
  brm_dat <- data.frame(X = X, Z = Z, Y1 = Y[,2], Y2 = Y[,3])
  # fit the model
  to_stan <- make_standata(mvbind(Y1, Y2) ~ X, data = brm_dat)
  m_mult <- stan("models/mult_model.stan",data = to_stan)
  # get the fitted residual correlation
  bias_m <- data.frame(Bias = median(extract(m_mult)$rescor - R[2,1]), Type = "model")
  bias_all <- rbind(bias_obs, bias_m)
  return(bias_all)
}

# replicate this 100 times
rep_exp7 <- replicate(100,sim_exp7(),simplify = FALSE)
rep_exp7 <- rbind.fill(rep_exp7)
rep_exp7$real_R <- -0.01
rep_exp7$fit_R <- with(rep_exp7, Bias + real_R)

ggplot(rep_exp7,aes(x=Type,y=fit_R)) +
  geom_jitter() +
  geom_hline(yintercept = rep_exp7$real_R, linetype = "dashed")

## eigth simulation, two response variables, two predictors, with differing effect on the responses and positive residual correlation, one predictors in the model
# the (residual) correlation matrix
R <- matrix(c(1, 0.6, 0.6, 1),ncol=2)
# the variable sd
stdev <- c(2, 1.5)
# get the residual covariance matrix
D <- cor2cov(R, stdev)
# the predictor
X <- runif(100, -2, 2)
Z <- runif(100, -2, 2)
# the linear predictor with different intercept for the two variables but same slope
mu <- cbind(1 + 2 * X + Z, 5 - 2 * X - Z)
# simulate the data packed in a function
sim_exp8 <- function(){
  Y <- adply(mu, 1, function(mu) rmvnorm(n = 1, mean = mu, sigma = D))
  # the observed correlation matrix
  bias_obs <- data.frame(Bias = cor(Y[,-1])[2,1] - R[2,1], Type = "observed")
  
  # format for brms
  brm_dat <- data.frame(X = X, Z = Z, Y1 = Y[,2], Y2 = Y[,3])
  # fit the model
  to_stan <- make_standata(mvbind(Y1, Y2) ~ X, data = brm_dat)
  m_mult <- stan("models/mult_model.stan",data = to_stan)
  # get the fitted residual correlation
  bias_m <- data.frame(Bias = median(extract(m_mult)$rescor - R[2,1]), Type = "model")
  bias_all <- rbind(bias_obs, bias_m)
  return(bias_all)
}

# replicate this 100 times
rep_exp8 <- replicate(100,sim_exp8(),simplify = FALSE)
rep_exp8 <- rbind.fill(rep_exp8)
rep_exp8$real_R <- 0.6
rep_exp8$fit_R <- with(rep_exp8, Bias + real_R)

ggplot(rep_exp8,aes(x=Type,y=fit_R)) +
  geom_jitter() +
  geom_hline(yintercept = rep_exp8$real_R, linetype = "dashed")

## nith simulation, two response variables, two predictors, with same effect on the responses and positive residual correlation, one predictors in the model
# the (residual) correlation matrix
R <- matrix(c(1, 0.6, 0.6, 1),ncol=2)
# the variable sd
stdev <- c(2, 1.5)
# get the residual covariance matrix
D <- cor2cov(R, stdev)
# the predictor
X <- runif(100, -2, 2)
Z <- runif(100, -2, 2)
# the linear predictor with different intercept for the two variables but same slope
mu <- cbind(1 + 2 * X + Z, 5 + 2 * X + Z)
# simulate the data packed in a function
sim_exp9 <- function(){
  Y <- adply(mu, 1, function(mu) rmvnorm(n = 1, mean = mu, sigma = D))
  # the observed correlation matrix
  bias_obs <- data.frame(Bias = cor(Y[,-1])[2,1] - R[2,1], Type = "observed")
  
  # format for brms
  brm_dat <- data.frame(X = X, Z = Z, Y1 = Y[,2], Y2 = Y[,3])
  # fit the model
  to_stan <- make_standata(mvbind(Y1, Y2) ~ X, data = brm_dat)
  m_mult <- stan("models/mult_model.stan",data = to_stan)
  # get the fitted residual correlation
  bias_m <- data.frame(Bias = median(extract(m_mult)$rescor - R[2,1]), Type = "model")
  bias_all <- rbind(bias_obs, bias_m)
  return(bias_all)
}

# replicate this 100 times
rep_exp9 <- replicate(100,sim_exp9(),simplify = FALSE)
rep_exp9 <- rbind.fill(rep_exp9)
rep_exp9$real_R <- 0.6
rep_exp9$fit_R <- with(rep_exp9, Bias + real_R)

ggplot(rep_exp9,aes(x=Type,y=fit_R)) +
  geom_jitter() +
  geom_hline(yintercept = rep_exp9$real_R, linetype = "dashed")
