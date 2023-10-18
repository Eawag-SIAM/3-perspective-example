## -------------------------------------------------------
##
## A toy example with a simple one-compartment toxicokinetic model
##
## October 18, 2023
## Sandrine Charles, Roman Ashauer, Andreas Scheidegger
##
## -------------------------------------------------------

library(ggplot2)
library(maxLik)            # maximum likelihood estimation
library(adaptMCMC)         # MCMC sampling algorithm

## Notes
## -----
## 
## - We use now the raw data from Ashauer 2010. No sure what 'samping_date' is.
##   (not the same as 'replicate' in Sandrins data set)
##
## - We use LOQ/2. This is ok for here but statistically not 100% clean.
##
## - TODO: compute independent priors based on Romans instuctions


pdf("plot.pdf", width=9)


## ---------------------------------
## 1) import data

## Data set on *Gammarus pulex* exposed to Malathion,  Ashauer (2010).
df <- read.table("data_Malathion_Ashauer_2010.csv", header = TRUE, sep = ",")
df$sampling_date <- as.factor(df$sampling_date)

ggplot(data = df, aes(x = time, y = Cinternal, col=sampling_date)) +
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
                   C_water,  # exposure concentration water
                   tc,       # duration of the accumulation phase [days]
                   tsim){    # simulation times

    k.u <- parameters[1]
    k.e <- parameters[2]
    
    ## accumlation phase
    tacc <- tsim[tsim <= tc]
    Cacc <- k.u*C_water*(1 - exp(-k.e*tacc))/k.e
    
    tdep <- tsim[tsim > tc]
    Cdep <- k.u*C_water*(exp(k.e*(tc - tdep)) - exp(-k.e*tdep))/k.e
    
    result <- data.frame(time = c(tacc, tdep),
                         conc = c(Cacc, Cdep))
    return(result)
}


## ---------------------------------
## 3) Perspective 1


## Paraemters form literature !!!CHECK!!!
ku.mean <- 92
ke.mean <-  0.815
parameters.P1 <- c(ku.mean, ke.mean)
C_water <- unique(df$Cwater)
tc <- 1 # duration of the accumulation phase
tsim <- seq(0, 7.5, length=111)


## run model
simu.mean <- bioacc(parameters.P1, C_water, tc, tsim)

## Plot model fit
ggplot(simu.mean, aes(time, conc)) +
    geom_line(col = "orange", linewidth = 1.5) +
    xlab("Time (hours)") +
    ylab("Internal predicted concentrations") +
    geom_vline(xintercept = 1, linetype="dashed") +    
    ## Add observed data
    geom_point(data = df, aes(time, Cinternal)) +
    ggtitle("Perspective 1")



## ---------------------------------
## 4) Perspective 2

## 4.1) --- log likelihood functions p(D|theta) assuming iid normal errors
##          theta_3 is the *log* of standard deviation of the error term
log.ll.normal <- function(theta, data, LOQ=6) {
    conc.pred <- bioacc(theta, unique(data$Cwater), tc = 1, tsim = data$time)$conc
    
    ## replace 0 observations with 1/2 LOQ
    obs <- ifelse(data$Cinternal==0, LOQ/2, data$Cinternal)
    
    sum(dnorm(obs, mean=conc.pred, sd=exp(theta[3]), log=TRUE))
}

log.ll.normal(c(k.u=90, k.e=0.8, log.sd=log(1)), df)

mle = maxLik(log.ll.normal,
             start=c(k.u=90, k.e=1, log.sd=log(1)),
             data=df,
             LOQ = LOQ)

summary(mle)


## run model
parameters.P2 <- mle$estimate
simu.mean <- bioacc(parameters.P2, C_water, tc, tsim)

## Plot model fit
ggplot(simu.mean, aes(time, conc)) +
    geom_line(col = "orange", linewidth = 1.5) +
    xlab("Time (hours)") +
    ylab("Internal predicted concentrations") +
    geom_vline(xintercept = 1, linetype="dashed") +    
    ## Add observed data
    geom_point(data = df, aes(time, Cinternal)) +
    ggtitle("MLE fit for normal errors")



## 4.2) --- log likelihood functions p(D|theta) assuming iid Gamma distributed
##          errors. The deterministic model determines the mode of the distribtion,
##          theta_3 the *log* of standard deviation of the error term.
log.ll.gamma <- function(theta, data, LOQ=6) {

    mode <- bioacc(theta, unique(data$Cwater), tc = 1, tsim = data$time)$conc
    sd <- exp(theta[3])

    ## calculate parameters of gamma distributions
    scale <- 0.5*(sqrt(4*sd^2+mode^2) - mode)
    shape <- mode/scale + 1

    ## replace 0 observations with 1/2 LOQ
    obs <- ifelse(data$Cinternal==0, LOQ/2, data$Cinternal)
    
    sum(dgamma(obs, shape=shape, scale=scale, log=TRUE))
}


mle = maxLik(log.ll.gamma,
             start=c(k.u=80, k.e=0.8, log.sd=log(1)),
             data=df,
             LOQ=LOQ,
             method="BFGS")
summary(mle)


## run model
parameters.P2 <- mle$estimate
simu.mean <- bioacc(parameters.P2, C_water, tc, tsim)

## Plot model fit
ggplot(simu.mean, aes(time, conc)) +
    geom_line(col = "orange", linewidth = 1.5) +
    xlab("Time (hours)") +
    ylab("Internal predicted concentrations") +
    geom_vline(xintercept = 1, linetype="dashed") +    
    ## Add observed data
    geom_point(data = df, aes(time, Cinternal)) +
    ggtitle("MLE fit for gamma errors")




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

log.post.normal <- function(theta, data, ...) {    
    log.ll.normal(theta, data, ...) + log.prior(theta)
}

log.post.gamma <- function(theta, data, ...) {    
    log.ll.gamma(theta, data, ...) + log.prior(theta)
}




## 5.3.1) --- sample from posterior with normal error

pp.normal <- MCMC(log.post.normal, n=15000,
                  init=c(k.u=90, k.e=1, log.sd=2),
           scale = c(100, 1, 1),
           acc.rate=0.2,
           data=df, LOQ = LOQ)
post.normal <- convert.to.coda(pp.normal)

plot(post.normal)                              # whole chain

## remove burn-in, and plot posterior distribution
plot(window(post.normal, start=5000))

pairs(pp.normal$samples, col=gray(0, alpha=0.01))


## 5.3.2) --- sample from posterior with Gamma error

pp.gamma <- MCMC(log.post.gamma, n=15000,
                 init=c(k.u=90, k.e=1, log.sd=2),
           scale = c(100, 1, 1),
           data=df, acc.rate=0.2)
post.gamma <- convert.to.coda(pp.gamma)

plot(post.gamma)                              # whole chain

## remove burn-in, and plot posterior distribution
plot(window(post.gamma, start=5000))

pairs(pp.gamma$samples, col=gray(0, alpha=0.01))



## 5.4.1) --- plot credibility intervals for normal errors

## sample from posterior
n <- 10000
X.para <- X.total <- matrix(NA, ncol=length(tsim), nrow=n)
for(i in 1:n){

    idx <- sample(1:nrow(pp.normal$samples),1)
    para <- pp.normal$samples[idx,]
    
    X.para[i,] <- bioacc(para, C_water, tc, tsim)$conc

    sd <- exp(para[3])
    X.total[i,] <- rnorm(length(tsim), mean=X.para[i,], sd=sd)

}

## compute quantiles
intervals.para <-  t(apply(X.para, 2,
                           function(x){quantile(x, probs=c(0.05, 0.5, 0.95))}))
intervals.para <- as.data.frame(intervals.para)
colnames(intervals.para) <- c("p05", "p50", "p95")
intervals.para$time <- tsim

intervals.total <-  t(apply(X.total, 2,
                            function(x){quantile(x, probs=c(0.05, 0.5, 0.95))}))
intervals.total <- as.data.frame(intervals.total)
colnames(intervals.total) <- c("p05", "p50", "p95")
intervals.total$time <- tsim


## Plot intervals
ggplot(intervals.para) +
    geom_line(aes(x=time, y=p50), col = "orange", linewidth = 1.5) +
    geom_ribbon(aes(x=time, y=p50, ymin=p05, ymax=p95), alpha = 0.2) +
    geom_ribbon(data=intervals.total, aes(x=time, y=p50, ymin=p05, ymax=p95), alpha = 0.2) +
    xlab("Time (hours)") +
    ylab("Internal concentration") +
    geom_vline(xintercept = 1, linetype="dashed") +    
    ## Add observed data
    geom_point(data = df, aes(time, Cinternal)) +
    ggtitle("90% parameter and total uncertainty interval with normal errors")



## 5.4.2) --- plot credibility intervals for gamma errors

## sample from posterior
n <- 10000
X.para <- X.total <- matrix(NA, ncol=length(tsim), nrow=n)
for(i in 1:n){

    idx <- sample(1:nrow(pp.gamma$samples),1)
    para <- pp.gamma$samples[idx,]
    
    X.para[i,] <- mode <- bioacc(para, C_water, tc, tsim)$conc
    sd <- exp(para[3])

    ## calculate parmeters of gamma distributions
    scale <- 0.5*(sqrt(4*sd^2+mode^2) - mode)
    shape <- mode/scale + 1
    
    X.total[i,] <- rgamma(length(tsim), scale=scale, shape=shape)

}

## compute quantiles
intervals.para <-  t(apply(X.para, 2,
                           function(x){quantile(x, probs=c(0.05, 0.5, 0.95))}))
intervals.para <- as.data.frame(intervals.para)
colnames(intervals.para) <- c("p05", "p50", "p95")
intervals.para$time <- tsim

intervals.total <-  t(apply(X.total, 2,
                            function(x){quantile(x, probs=c(0.05, 0.5, 0.95))}))
intervals.total <- as.data.frame(intervals.total)
colnames(intervals.total) <- c("p05", "p50", "p95")
intervals.total$time <- tsim


## Plot credibility intervals
ggplot(intervals.para) +
    geom_line(aes(x=time, y=p50), col = "orange", linewidth = 1.5) +
    geom_ribbon(aes(x=time, y=p50, ymin=p05, ymax=p95), alpha = 0.2) +
    geom_ribbon(data=intervals.total, aes(x=time, y=p50, ymin=p05, ymax=p95), alpha = 0.2) +
    xlab("Time (hours)") +
    ylab("Internal concentration") +
    geom_vline(xintercept = 1, linetype="dashed") +    
    ## Add observed data
    geom_point(data = df, aes(time, Cinternal)) +
    ggtitle("90% parameter and total uncertainty interval with gamma errors")


dev.off()
