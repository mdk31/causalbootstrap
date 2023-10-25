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
ns <- c(5e3, 10e3, 20e3)
z_a <- qnorm(0.975)
treat_eff <- 0.5

image_path <- '~/Documents/HW/Research/CI/causalbootstrap/images'
dat_path <- '~/Documents/HW/Research/CI/causalbootstrap/data'

setwd('~/Documents/HW/Research/CI/causalbootstrap/code')

blb_params <- expand.grid(B = c(100, 500),
                          subsets = c(2, 4, 10),
                          method = c('svmLinear2', 'ps', 'cbps'),
                          n_size = ns,
                          sigma = c(1),
                          gamma = c(0.7, 0.75, 0.8),
                          te = c(treat_eff),
                          dgp_func = c('kangschafer3_dep', 'kangschafer3', 'kangschafer3_mult',
                                       'kangschafer3_het', 'kangschafer3_mis'),
                          stringsAsFactors = FALSE)

combos <- c('n', 'sigma', 'te', 'B', 'subsets', 'method', 'gamma', 'dgp_func', 'boot_type')

source('causal_funcs.R')

if(!file.exists(file.path(dat_path, 'cblb_sims.rds'))){
  cblb <- lapply(seq_len(nrow(blb_params)), function(didx){
    dg_row <- blb_params[didx, ]
    print(dg_row)
    if(dg_row$dgp_func == 'kangschafer3_mult'){
      form <- Tr ~ 1 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12
      overlap <- 0.3
    } else{
      form <- Tr ~ 1 + X1 + X2
      overlap <- 0.5
    }
    replicates <- pblapply(seq_len(reps), function(repnum){
      set.seed(repnum)
      dat <- do.call(dg_row$dgp_func, args = list(dg_row$n_size, te = dg_row$te, 
                                                  sigma = dg_row$sigma, beta_overlap = overlap))
      # Weighted BLB
      part <- weight_blb(data = dat, 
                         R = dg_row$B, 
                         pi_formula = as.formula(form), 
                         Y = 'y', 
                         W = 'Tr',
                         b = round(nrow(dat)^dg_row$gamma),
                         subsets = dg_row$subsets,
                         method = dg_row$method,
                         type = 'obs',
                         prop_grid = data.frame(cost = 0.01))
      
      part <- purrr::transpose(part)
      theta_reps <- part$theta_reps
      
      perc_ci <- lapply(theta_reps, function(boot_reps){
        boot:::perc.ci(t = boot_reps)[, 4:5]
      })
      perc_ci <- Reduce(`+`, perc_ci)/dg_row$subsets
      # NORMAL
      norm_ci <- purrr::map2(theta_reps, part$theta0, ~{
        boot:::norm.ci(var.t0 = NULL, t0 = .y, t = .x)[, 2:3]
      })
      norm_ci <- Reduce(`+`, norm_ci)/dg_row$subsets
      
      estim <- sapply(theta_reps, mean)
      estim <- mean(estim)
      se <- sapply(theta_reps, sd)
      se <- mean(se)
      
      data.table(estim = estim,
                 se = se,
                 lower = c(perc_ci[1], norm_ci[1]),
                 upper = c(perc_ci[2], norm_ci[2]),
                 B = dg_row$B,
                 subsets = dg_row$subsets,
                 method = dg_row$method,
                 truth = unlist(part$truth)[1],
                 boot_type = c('percentile', 'normal'))
      
    }, cl = 4)
    out <- rbindlist(replicates)
    out[, `:=`(n = dg_row$n_size,
               sigma = dg_row$sigma,
               gamma = dg_row$gamma,
               te = dg_row$te,
               dgp_func = dg_row$dgp_func)]
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(dat_path, 'cblb_sims.rds'))
} else{
  cblb <- readRDS(file.path(dat_path, 'cblb_sims.rds'))
}

# BIAS (LOOK AT BASELINE WITH GAMMA = 0.8)
ggdat <- copy(cblb)
ggdat <- ggdat[boot_type == 'percentile' & dgp_func == 'kangschafer3' & gamma == 0.8]
ggdat[, `:=`(n = factor(n))]
ggdat[, `:=`(mn = mean(estim)), by = c(combos)]
ggdat[, `:=`(ystrt = as.integer(n) - 0.25,
             yend = as.integer(n) + 0.25,
             B = factor(B),
             method = car::recode(method, "'ps'='PS';'cbps'='CBPS';'svmLinear2'='SVM'"))]
p <- ggplot(ggdat, aes(x = estim, y = n, color = B), alpha = 0.2) +
  geom_point() +
  geom_jitter(width = 0, height = 0.25) +
  facet_grid( method ~ subsets, labeller = label_both) +
  geom_vline(aes(xintercept = te), color = 'yellow', size = 0.5, linetype = 'dashed') +
  geom_segment(aes(x = mn, y = ystrt, xend = mn, yend = yend),
               color = 'red', linetype = 'dashed', size = 0.5) +
  ylab('Sample size') +
  xlab('ATE Estimate') +
  xlim(treat_eff + 0.75*c(-1, 1)) +
  theme_bw()
print(p)
ggsave(file.path(image_path, 'blb_bias.pdf'), height = 9, width = 7)

# COVERAGE----
cblb[, `:=`(dgp_func = car::recode(dgp_func, "'kangschafer3_dep'='Dependent Variables';'kangschafer3'='Baseline';
                                   'kangschafer3_mult'='Multiple Predictors';'kangschafer3_het'='Heterogeneity';
                                   'kangschafer3_mis'='Misspecification'"))]
zip <- copy(cblb)
zip <- zip[B == 100]
zip[, `:=`(cent = abs(estim - te)/se)]
zip[, `:=`(rank = as.integer(cut(cent, quantile(cent, probs = seq(0, 1, by = 0.01)), include.lowest = TRUE))), by = combos]

zip[, `:=`(n = format(n, scientific = FALSE, big.mark = ','),
           covered = fifelse(lower <= te & upper >= te, 'Coverer', 'Non-coverer'))]
zip_labels <- zip[, .(perc_cover = round(mean(covered == 'Coverer'), 3)),
                  by = combos]

dgp_fix <- unique(cblb$dgp_func)
sub_fix <- unique(cblb$subsets)
gamma_fix <- unique(cblb$gamma)
ci_fix <- unique(cblb$boot_type)

ggdat <- copy(zip)
for(i in sub_fix){
  for(j in dgp_fix){
    for(k in gamma_fix){
      for(l in ci_fix){
        nm <- paste0('blb_zip_sub', i, '_', paste(j, collapse = '_'), '_gamma_', k, '_ci_', l, '.pdf')
        title <- bquote(paste(s == .(i), ' and ', gamma == .(k), ' and ', .(j)))
        ggsub <- ggdat[subsets == i & dgp_func == j & gamma == k & boot_type == l]
        label_sub <- zip_labels[subsets == i & dgp_func == j & gamma == k & boot_type == l]
        p <- ggplot(ggsub, aes(y = rank)) +
          geom_segment(aes(x = lower, y = rank, xend = upper, yend = rank, color = covered)) +
          facet_grid(method ~ n, labeller = label_both) +
          geom_vline(aes(xintercept = te), color = 'yellow', size = 0.5, linetype = 'dashed') +
          ylab('Fractional Centile of |z|') +
          xlab('95% Confidence Intervals') +
          theme_bw() +
          scale_y_continuous(breaks = c(5, 50, 95)) +
          scale_color_discrete(name = "Coverage") +
          geom_text(x = 0.65, y = 50, aes(label = perc_cover), data = label_sub, size = 4) +
          ggtitle(title)
        
        print(p)
        ggsave(file.path(image_path, nm), height = 9, width = 7)
      }

    }
  }
}



