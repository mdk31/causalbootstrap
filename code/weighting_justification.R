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

n <- 5000
te <- 0.5
r <- 100
b <- round(n^0.8)
idx <- seq_len(n)
part_idx <- seq_len(b)

image_path <- '~/Documents/HW/Research/CI/causalbootstrap/images'
dat_path <- '~/Documents/HW/Research/CI/causalbootstrap/data'

setwd('~/Documents/HW/Research/CI/causalbootstrap/code')

source('causal_funcs.R')

dat <- kangschafer3(n = 5000, te = te, sigma = 1)
part_dat <- dat[sample(idx, size = b, replace = FALSE)]

wts <- make_weights(Tr ~ 1 + X1 + X2, method = 'ps', data = part_dat, normed = TRUE)
part_dat$wts <- wts

boot_reps_n <- replicate(r, {
  boot_idx <- sample(part_idx, size = n, replace = TRUE)
  boot_dat <- part_dat[boot_idx]
  ipw <- (boot_dat$Tr)*boot_dat$y*boot_dat$wts - (1 - boot_dat$Tr)*boot_dat$y*boot_dat$wts
  sum(ipw)
})

boot_reps_b <- replicate(r, {
  boot_idx <- sample(part_idx, size = b, replace = TRUE)
  boot_dat <- part_dat[boot_idx]
  ipw <- (boot_dat$Tr)*boot_dat$y*boot_dat$wts - (1 - boot_dat$Tr)*boot_dat$y*boot_dat$wts
  sum(ipw)
})

reps <- data.table(`Bootstrap Size` = rep(c('b', 'n'), each = r),
           estim = c(boot_reps_b, boot_reps_n),
           truth = 0.5)

p <- ggplot(reps, aes(x = estim)) +
  geom_histogram() +
  facet_wrap(~ `Bootstrap Size`) +
  geom_vline(aes(xintercept = te), color = 'yellow', size = 0.5, linetype = 'dashed') +
  xlab('ATE Estimates') +
  theme_bw()

print(p)
ggsave(file.path(image_path, 'weight_justification.pdf'), height = 9, width = 7)


