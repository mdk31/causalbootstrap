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

n <- 10000
te <- 0.5
r <- 100
idx <- seq_len(n)

image_path <- '~/Documents/HW/Research/CI/causalbootstrap/images'
dat_path <- '~/Documents/HW/Research/CI/causalbootstrap/data'

setwd('~/Documents/HW/Research/CI/causalbootstrap/code')

source('causal_funcs.R')


out <- replicate(500, {
  dat <- kangschafer3(n = n, te = te, sigma = 1)
  wts <- WeightIt::weightit(formula = Tr ~ 1 + X1 + X2,
                            data = dat,
                            method = 'ps',
                            estimand = 'ATE',
                            over = FALSE)
  dat$wts <- wts$weights
  boot_reps_n <- replicate(r, {
    boot_idx <- sample(idx, size = n, replace = TRUE)
    boot_dat <- dat[boot_idx]
    ipw <- (boot_dat$Tr)*boot_dat$y*boot_dat$wts - (1 - boot_dat$Tr)*boot_dat$y*boot_dat$wts
    mean(ipw)
  })

  boot_ci <- boot:::perc.ci(t = boot_reps_n)[, 4:5]
  data.table::between(te, boot_ci[1], boot_ci[2])
})
cat('Logistic regression nonparametric bootstrap: ', mean(out))

out <- replicate(500, {
  dat <- kangschafer3(n = n, te = te, sigma = 1)
  wts <- make_weights(Tr ~ 1 + X1 + X2, data = dat, method = 'cbps', normed = FALSE)
  dat$wts <- wts
  boot_reps_n <- replicate(r, {
    boot_idx <- sample(idx, size = n, replace = TRUE)
    boot_dat <- dat[boot_idx] 
    ipw <- (boot_dat$Tr)*boot_dat$y*boot_dat$wts - (1 - boot_dat$Tr)*boot_dat$y*boot_dat$wts
    mean(ipw)
  })
  
  boot_ci <- boot:::perc.ci(t = boot_reps_n)[, 4:5]
  data.table::between(te, boot_ci[1], boot_ci[2])
})
cat('CBPS nonparametric bootstrap: ', mean(out))



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


