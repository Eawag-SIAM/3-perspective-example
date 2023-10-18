
## -------------------------------------------------------
##
## A toy example with a simple one-compartment toxicokinetic model
##
## August 25, 2023
## Sandrine Charles, Roman Ashauer, Andreas Scheidegger
##
## -------------------------------------------------------

library(ggplot2)
library(maxLik)            # maximum likelihood estimation
library(adaptMCMC)         # MCMC sampling algorithm

## TODO/Questions
##
## - how can the final concentration be lower than the initial concentration?
##    -> is time 0 really 0? Check with raw data.
## - expw == exposer concentration in water --> rename
## - add confidence and prediction intervals
##
## perspective 1:
## - new values from Roman
##
## perspective 2:
## - how to deal with 0 concentrations for gamma error model? 
##    -> use 12 LOQ from Paper
##
## perspective 3
## - use prior from persp 1



## ---------------------------------
## 1) import data

## Data set on *Gammarus pulex* exposed to Malathion,  Ashauer (2010).
df <- read.table("data.txt", header = TRUE, sep = "")
df$replicate <- as.factor(df$replicate)

ggplot(data = df, aes(x = time, y = conc, col=replicate)) +
    geom_line() +
    geom_point(size=5) +
    xlab("Time (days)") +
    ylab("Internal measured concentration (picomol/g wet weight)") +
    geom_vline(xintercept = 1, linetype="dashed")


## Limit of detection nmol/kg(wet weight)
## from Table 5 in SI of Ashaue et al (2010)
LOQ <- 6


## ---------------------------------
## 2) define model

## 2.1) --- the deterministic model with two parameters

bioacc <- function(parameters,
                   expw,  # ? exposure concentration water ?
                   tc,    # duration of the accumulation phase [h]
                   tsim){ # simulation times

    k.u <- parameters[1]
    k.e <- parameters[2]
    
    ## accumlation phase
    tacc <- tsim[tsim <= tc]
    Cacc <- k.u*expw*(1 - exp(-k.e*tacc))/k.e
    
    tdep <- tsim[tsim > tc]
    Cdep <- k.u*expw*(exp(k.e*(tc - tdep)) - exp(-k.e*tdep))/k.e
    
    result <- data.frame(time = c(tacc, tdep),
                         conc = c(Cacc, Cdep))
    return(result)
}


## ---------------------------------
## 3) Perspective 1


## parameters values from Ashauer et al. 2010

para.ashauer <- read.table("parameters_Ashauer.csv",
                           header = TRUE, sep=";")

## Assign input parameter values
## Use the mean of Marquardt and MCMC fits
ku.mean <- mean(para.ashauer$ku_mean) 
ke.mean <- mean(para.ashauer$ke_mean)

parameters.P1 <- c(ku.mean, ke.mean)
expw <- unique(df$expw)
tc <- 1 # duration of the accumulation phase
tsim <- seq(0, 7.5, length=111)


## run model
simu.mean <- bioacc(parameters.P1, expw, tc, tsim)

## Plot model fit
ggplot(simu.mean, aes(time, conc)) +
    geom_line(col = "orange", linewidth = 1.5) +
    xlab("Time (hours)") +
    ylab("Internal predicted concentrations") +
    geom_vline(xintercept = 1, linetype="dashed") +    
    ## Add observed data
    geom_point(data = df, aes(time, conc))



## ---------------------------------
## 4) Perspective 2

## 4.1) --- log likelihood functions p(D|theta) assuming iid normal errors
##          theta_3 is the *log* of standard deviation of the error term
log.ll.normal <- function(theta, data) {
    conc.pred <- bioacc(theta, unique(data$expw), tc = 1, tsim = data$time)$conc
    sum(dnorm(data$conc, mean=conc.pred, sd=exp(theta[3]), log=TRUE))
}


mle = maxLik(log.ll.normal, start=c(k.u=90, k.e=0.8, log.sd=log(1)), data=df)

summary(mle)


## run model
parameters.P2 <- mle$estimate
simu.mean <- bioacc(parameters.P2, expw, tc, tsim)

## Plot model fit
ggplot(simu.mean, aes(time, conc)) +
    geom_line(col = "orange", linewidth = 1.5) +
    xlab("Time (hours)") +
    ylab("Internal predicted concentrations") +
    geom_vline(xintercept = 1, linetype="dashed") +    
    ## Add observed data
    geom_point(data = df, aes(time, conc))



## 4.2) --- log likelihood functions p(D|theta) assuming iid Gamma distributed
##          errors. The deterministic model determines the mode of the distribtion,
##          theta_3 the *log* of standard deviation of the error term.
log.ll.gamma <- function(theta, data) {

    mode <- bioacc(theta, unique(data$expw), tc = 1, tsim = data$time)$conc
    sd <- exp(theta[3])

    ## calculate parameters of gamma distributions
    scale <- 0.5*(sqrt(4*sd^2+mode^2) - mode)
    shape <- mode/scale + 1

    ## replace 0 observations with LOQ
    obs <- ifelse(data$conc==0, 0.1, data$conc) 
    sum(dgamma(obs, shape=shape, scale=scale, log=TRUE))
}


mle = maxLik(log.ll.gamma,
             start=c(k.u=80, k.e=0.8, log.sd=log(1)),
             data=df,
             method="BFGS")
summary(mle)


## run model
parameters.P2 <- mle$estimate
simu.mean <- bioacc(parameters.P2, expw, tc, tsim)

## Plot model fit
ggplot(simu.mean, aes(time, conc)) +
    geom_line(col = "orange", linewidth = 1.5) +
    xlab("Time (hours)") +
    ylab("Internal predicted concentrations") +
    geom_vline(xintercept = 1, linetype="dashed") +    
    ## Add observed data
    geom_point(data = df, aes(time, conc))




## ---------------------------------
## 5) Perspective 3


## 5.1) --- define log prior !!!CHECK!!!
log.prior <- function(theta) {
    dnorm(theta[1], 90, 10, log=TRUE) +
        dunif(theta[2], 0, 2, log=TRUE) +
        dnorm(theta[3], 2, 2, log=TRUE) 
}


## 5.2) --- define function proportional to posterior distribution
##          p(theta | D) = p(D|theta) * p(theta) / p(D)

log.post.normal <- function(theta, data) {    
    log.ll.normal(theta, data) + log.prior(theta)
}

log.post.gamma <- function(theta, data) {    
    log.ll.gamma(theta, data) + log.prior(theta)
}




## 5.3.1) --- sample from posterior with normal error

pp <- MCMC(log.post.normal, n=15000,
           init=c(k.u=90, k.e=1, log.sd=2),
           scale = c(100, 1, 1),
           data=df, acc.rate=0.2)
post <- convert.to.coda(pp)

plot(post)                              # whole chain

## remove burn-in, and plot posterior distribution
plot(window(post, start=5000))


## 5.3.2) --- sample from posterior with Gamma error

pp <- MCMC(log.post.gamma, n=15000,
           init=c(k.u=90, k.e=1, log.sd=2),
           scale = c(100, 1, 1),
           data=df, acc.rate=0.2)
post <- convert.to.coda(pp)

plot(post)                              # whole chain

## remove burn-in, and plot posterior distribution
plot(window(post, start=5000))

pairs(pp$samples, col=gray(0, alpha=0.01))



## 5.4.1) --- plot credibility intervals for normal errors

n <- 10000
X <- matrix(NA, ncol=length(tsim), nrow=n)
for(i in 1:n){

    idx <- sample(1:nrow(pp$samples),1)
    para <- pp$samples[idx,]
    
    mean <- bioacc(para, expw, tc, tsim)$conc
    sd <- exp(para[3])

    
    X[i,] <- rnorm(length(tsim), mean=mean, sd=sd)

}

intervals <- as.data.frame(cbind(tsim,
                                 t(apply(X, 2,
                                         function(x){quantile(x, probs=c(0.05, 0.5, 0.95))}))))
colnames(intervals) <- c("time", "p05", "p50", "p95")


## Plot credibility intervals
ggplot(intervals) +
    geom_line(aes(x=time, y=p50), col = "orange", linewidth = 1.5) +
    geom_ribbon(aes(x=time, y=p50, ymin = p05, ymax = p95), alpha = 0.2) +
    xlab("Time (hours)") +
    ylab("Internal predicted concentrations") +
    geom_vline(xintercept = 1, linetype="dashed") +    
    ## Add observed data
    geom_point(data = df, aes(time, conc)) +
    ggtitle("90% credibility interval with normal errors")



## 5.4.2) --- plot credibility intervals for gamma errors

n <- 10000
X <- matrix(NA, ncol=length(tsim), nrow=n)
for(i in 1:n){

    idx <- sample(1:nrow(pp$samples),1)
    para <- pp$samples[idx,]
    
    mode <- bioacc(para, expw, tc, tsim)$conc
    sd <- exp(para[3])

    ## calculate parmeters of gamma distributions
    scale <- 0.5*(sqrt(4*sd^2+mode^2) - mode)
    shape <- mode/scale + 1
    
    X[i,] <- rgamma(length(tsim), scale=scale, shape=shape)

}

intervals <- as.data.frame(cbind(tsim,
                                 t(apply(X, 2,
                                        function(x){quantile(x, probs=c(0.05, 0.5, 0.95))}))))
colnames(intervals) <- c("time", "p05", "p50", "p95")


## Plot credibility intervals
ggplot(intervals) +
    geom_line(aes(x=time, y=p50), col = "orange", linewidth = 1.5) +
    geom_ribbon(aes(x=time, y=p50, ymin = p05, ymax = p95), alpha = 0.2) +
    xlab("Time (hours)") +
    ylab("Internal predicted concentrations") +
    geom_vline(xintercept = 1, linetype="dashed") +    
    ## Add observed data
    geom_point(data = df, aes(time, conc)) +
    ggtitle("90% credibility interval with gamma errors")

