## -------------------------------------------------------
##
## Derive priors for k.e and k.u bases on 
##
## October 18, 2023 -- Andreas Scheidegger
## andreas.scheidegger@eawag.ch
## -------------------------------------------------------

## CHECK:
## - good source of K.OW?
## - what is the base of reported log K.OW?
## - how should we guess the uncertainty?

## Based on:
## Hendriks, A. J., van der Linde, A., Cornelissen, G., & Sijm,
## D. T. H. M. (2001). The power of size. 1. Rate constants and
## equilibrium ratios for accumulation of organic substances related to
## octanol-water partition ratio and species weight. Environmental
## Toxicology and Chemistry, 20(7),
## 1399–1420. https://doi.org/10.1002/etc.5620200703

## -----------
## Constants, from Roman

w <- 2.25e-5        # organsim weight [kg]
kappa <- 0.25       # rate exponent [-]
rho.H2O <- 2.8e-3   # water layer diffusion rate [d kg^-kappa]
rho.CH2 <- 68       # lipid layer diffusion rate [d kg^-kappa]
p.CH2 <- 0.013      # Lipid fraction of organism [-]
gamma.0 <- 200      # Water absorption–excretion coefficient [gf^kappa/d]

## source: http://npic.orst.edu/factsheets/archive/malatech.html
K.OW <- 10^(2.75)   # Octanol–water partition ratio of Malathion [-]


## -----------
## Uptake rate k.u. Bases on Eq(5) of Hendriks et al.

k.u <- w^(-kappa) / (rho.H2O + rho.CH2/K.OW + 1/gamma.0)
k.u

## r2 on log scale: 0.39  (Fig 2, left)
##  -> uncertainty about +/- one order of magnitude

## -----------
## Elimination rate k.e. Based on Eq(8) of Hendriks et al.

k.e <- 1/(rho.CH2 * (K.OW - 1) + 1) * w^(-kappa)/(rho.H2O + rho.CH2/K.OW + 1/gamma.0)
k.e

## r2 on log scale: 0.7  (Fig 2, right) ???
## -> uncertainty about +/- one order of magnitude


## ----------
## Define distribution
## We need to have a guess about the uncertainty


## -- Idea I
## Use log normal with a mean of the above computation
## and relative standard deviation of one order of magnitude
## --> give us a extremely skewd distribution

mean <- k.u
sd <- 10*mean

meanlog <- log(mean) - 0.5*log(1 + (sd/mean)^2)
sdlog <- sqrt(log(1 + sd^2/(mean^2)))

x <- rlnorm(1000000, meanlog, sdlog)
mean(x)
sd(x)
hist(x)
hist(log10(x))
abline(v=log10(k.u), col=2)

plot(function(x) dlnorm(x, meanlog, sdlog), xlim=c(0, 500), n=301)
abline(v=k.u, col=2)


## -- Idea II
## Assume k.u is unbiased normal in the log space.
## -> a bit less skewed. median (k.u) = 10^log(mean k.u) 

meanlog = log10(k.u) / log10(exp(1))
sdlog = log10(10) / log10(exp(1)) # one order of magnitude (10 base)


x <- rlnorm(1000000, meanlog, sdlog)

median(x)
sd(x)
hist(x)
hist(log10(x))
abline(v=log10(k.u), col=2)

plot(function(x) dlnorm(x, meanlog, sdlog), xlim=c(0, 500), n=301)
abline(v=k.u, col=2)

