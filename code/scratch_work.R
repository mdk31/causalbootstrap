library(purrr)
library(data.table)
library(mice)
library(ggplot2)
library(dplyr)
library(kableExtra)
library(microbenchmark)
library(pbapply)
library(WeightIt)
set.seed(123)

# Microbenchmark replications
ntimes <- 100L
reps <- 500L
# ns <- c(5e3, 10e3, 20e3)
z_a <- qnorm(0.975)

te <- 0.5
sigma <- 1
n <- 5000
r <- 500

image_path <- '~/Documents/HW/Research/CI/causalbootstrap/images'
dat_path <- '~/Documents/HW/Research/CI/causalbootstrap/data'

setwd('~/Documents/HW/Research/CI/causalbootstrap/code')
source('causal_funcs.R')

idx <- seq_len(n)
W <- 'Tr'
Y <- 'y'

out <- replicate(reps, {
  dat <- kangschafer3_mult(n = n, sigma = sigma, te = te, preds = 2)
  n1 <- sum(dat$Tr == 1)
  n0 <- sum(dat$Tr == 0)
  
  wts <- make_weights(formula = Tr ~ 1 + V3 + V4,
                     data = dat,
                     method = 'ps',
                     normed = TRUE)

  dat$wts1 <- wts

  boot_out <- replicate(r, {
    boot_idx <- sample(idx, size = n, replace = TRUE)
    boot_dat <- dat[boot_idx]
    with(boot_dat, sum(Tr*y*wts1 - (1-Tr)*y*wts1))
  })
  cis <- quantile(boot_out, probs = c(0.025, 0.975))

  # ws <- split(wts, dat[[W]])
  # ys <- split(dat[[Y]], dat[[W]])
  # w0 <- ws[['0']]
  # w1 <- ws[['1']]
  # 
  # M1 <- rmultinom(r, n1, prob = w1)
  # y1 <- colSums(ys[['1']]*M1)/n1
  # M0 <- rmultinom(r, n0, prob = w0)
  # y0 <- colSums(ys[['0']]*M0)/n0
  # 
  # cis <- quantile(y1 - y0, probs = c(0.025, 0.975))
  data.table::between(te, cis[1], cis[2])
  
})
mean(out)
