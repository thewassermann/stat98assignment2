# load libraries
library(data.table)
library(functional)
library(lmtest)

# set seed for reproducability
set.seed(98)

# This is the script for Simulation Project for Group D, STAT 98 Spring Semester 2017
# Group Members: Nathan Goldberg, Alyssa Mehta and Gil Wassermann

total.func <- function(n.sims, n, line.params, x.gen, error.func, sig.level){
  # error.func : function, takes x as an input and returns
  # a random error value based on a distribution
  # n.sims: int, number of iterations 
  # n : int, number of data points to generate per simulation
  # x.gen : function,  generates the values of x for the data using a distribution
  # line.params : vector of doubles, parameters governing the linear line on which the
  # simulation is based
  # sig.level : double,  for the CI for coverage probability and BP test
  
  # actual parameters
  BETA0 <- line.params[1]
  BETA1 <- line.params[2]
  
  # create a data.table for results
  res.table <- data.table(
    b0.mle=rep(NA, n.sims),
    b0.cover.bool=rep(NA, n.sims),
    b0.ci.lwr=rep(NA, n.sims),
    b0.ci.upr=rep(NA, n.sims),
    b1.mle=rep(NA, n.sims),
    b1.cover.bool=rep(NA, n.sims),
    b1.ci.lwr=rep(NA, n.sims),
    b1.ci.upr=rep(NA, n.sims),
    s2.mle=rep(NA, n.sims),
    bp.test=rep(NA, n.sims)
  )  
  
  # loop through number of simulations
  for (i in 1:n.sims){
    
    # generate xs and ys 
    xy.DT <- sim.func(n, x.gen, error.func, line.params)
    
    # fit an OLS regression
    fit <- lm(xy.DT$y ~ xy.DT$x)
    
    # summary of data
    OLS.summary <- summary(fit)
    
    # breusch-pagan test
    bp.res <- bptest(fit)
    
    # confidence intervals
    param.confint <- confint(fit, level=sig.level)
    
    # extract and fill estimator and ci data
    res.table$b0.mle[i] <- OLS.summary$coefficients[1,1]
    res.table$b1.mle[i] <- OLS.summary$coefficients[2,1]
    res.table$s2.mle[i] <- OLS.summary$sigma**2
    
    res.table$b0.ci.lwr[i] <-param.confint[1,1]
    res.table$b0.ci.upr[i] <-param.confint[1,2]
    res.table$b1.ci.lwr[i] <-param.confint[2,1]
    res.table$b1.ci.upr[i] <-param.confint[2,2]
    
    # fill boolean data
    if ((BETA0 >= res.table$b0.ci.lwr[i]) & (BETA0 <= res.table$b0.ci.upr[i])){
      res.table$b0.cover.bool[i] <- TRUE
    } else {
      res.table$b0.cover.bool[i] <- FALSE 
    }
    
    if ((BETA1 >= res.table$b1.ci.lwr[i]) & (BETA1 <= res.table$b1.ci.upr[i])){
      res.table$b1.cover.bool[i] <- TRUE  
    } else {
      res.table$b1.cover.bool[i] <- FALSE  
    }
    
    # add hyp test results
    if( bp.res$p.value < 1-sig.level){
      res.table$bp.test[i] <- 'reject'
    } else {
      res.table$bp.test[i] <- 'accept'
    }
    
  }
  return(res.table)
}

sim.func <- function(n, x.gen, error.func, line.params){
  # function for a single iteration of the simulation
  
  # create data.table to store simulated data
  sim.table <- data.table(
    x=rep(NA, n),
    y=rep(NA, n)
  )
  
  # extract values from line.params
  beta.0 <- line.params[1]
  beta.1 <- line.params[2]
  
  for (i in 1:n){
    
    # pick a value of X from a x.gen distribution
    x <- x.gen(n=1)
    err <- error.func(x)
    
    # get y value
    y <- beta.0 + (beta.1 * x) + err
    
    # add value to sim.table
    sim.table$x[i] <- x
    sim.table$y[i] <- y
  }
  
  # return `n` simulated values of
  # x and y
  return(sim.table)
}

save.table <- function(tbl, fname){
  # function to save results to disc
  # tbl : data.table
  # fname : string for file name
  write.csv(tbl, file = fname)
}

###############
# Test Script #
###############
# x.gen.test <- Curry(runif)
# error.func.test <- function(x) return(Curry(rnorm, n=1, sd=1)())
# res <- total.func(100, 10, c(0,2), x.gen.test, error.func.test, 0.95)


##############
# Simulation #
##############

# params for simulation
N.SIMS <- 10000
N <- 100

# X generation function #
x.gen <- Curry(runif)

# ERROR FUNCTIONS #

# ~ norm(0, var=X_i)
error.norm.x <- function(x) return(Curry(rnorm, n=1, sd=sqrt(x))())

# ~ norm(0, sigma2) where sigma2 ~ Gamma(X_i, 1)
error.norm.gamma1 <- function(x){
  sigma <- sqrt(rgamma(1, x, 1))
  return(Curry(rnorm, n=1, sd=sigma)())
}

# ~ norm(0, sigma2) where sigma2 ~ Gamma(1, X_i)
error.norm.gamma2 <- function(x){
  sigma <- sqrt(rgamma(1, 1, x))
  return(Curry(rnorm, n=1, sd=sigma)())
}

# ~unif(-x/2, x/2)
error.unif.x <- function(x) return(Curry(runif, n=1, min=(-x/2), max=(x/2))())

# Create vector of functions for looping
func.vect <- c(error.norm.x, error.norm.gamma1, error.norm.gamma2, error.unif.x)

# loop for full simulation
for (i in 1:length(func.vect)){
  
  # name for file
  fname <- paste('sim_', i, sep='')
  fname <- paste(fname, '.csv', sep='')
  
  # run simulation
  sim.res <- total.func(N.SIMS, N, c(0,2), x.gen, func.vect[[i]], 0.95)
  
  # print to csv
  save.table(sim.res, fname)
}
