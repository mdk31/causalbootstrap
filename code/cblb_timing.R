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
ns <- c(5e3, 10e3, 20e3)
z_a <- qnorm(0.975)
treat_eff <- 0.5

image_path <- '~/Documents/HW/Research/CI/causalbootstrap/images'
dat_path <- '~/Documents/HW/Research/CI/causalbootstrap/data'

setwd('~/Documents/HW/Research/CI/causalbootstrap/code')

blb_params <- expand.grid(B = c(100),
                          subsets = c(2, 4, 10),
                          method = c('svmLinear2', 'ps', 'cbps'),
                          n_size = ns,
                          sigma = c(1),
                          gamma = c(0.7, 0.75, 0.8),
                          te = c(treat_eff),
                          dgp_func = c('kangschafer3'),
                          stringsAsFactors = FALSE)
tmp_params <- expand.grid(B = c(100),
                          subsets = c(1),
                          method = c('svmLinear2', 'ps', 'cbps'),
                          n_size = ns,
                          sigma = c(1),
                          gamma = c(1),
                          te = c(treat_eff),
                          dgp_func = c('kangschafer3'),
                          stringsAsFactors = FALSE)
blb_params <- rbind(blb_params, tmp_params)

combos <- c('n', 'sigma', 'te', 'B', 'subsets', 'method', 'gamma')

source('causal_funcs.R')

if(!file.exists(file.path(dat_path, 'cblb_timing.rds'))){
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
    replicates <- lapply(seq_len(ntimes), function(repnum){
      set.seed(repnum)
      dat <- do.call(dg_row$dgp_func, args = list(dg_row$n_size, te = dg_row$te, 
                                                  sigma = dg_row$sigma, beta_overlap = overlap))
      # Weighted BLB
      time_out <- system.time({
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
        
        estim <- sapply(theta_reps, mean)
        estim <- mean(estim)
        se <- sapply(theta_reps, sd)
        se <- mean(se)
      })

      data.table(time_out = time_out['elapsed'])
      
    })
    out <- rbindlist(replicates)
    out[, `:=`(n = dg_row$n_size,
               B = dg_row$B,
               subsets = dg_row$subsets,
               sigma = dg_row$sigma,
               gamma = dg_row$gamma,
               method = dg_row$method,
               te = dg_row$te)]
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(dat_path, 'cblb_timing.rds'))
} else{
  cblb <- readRDS(file.path(dat_path, 'cblb_timing.rds'))
}

cblb[, `:=`(subsets = as.character(subsets))]
cblb[, `:=`(subsets = fifelse(subsets == 1, 'Full Sample', subsets))]
cblb[, `:=`(subsets = factor(subsets, levels = c( '2', '4', '10', 'Full Sample')))]
ggdat <- copy(cblb)
ggdat <- ggdat[gamma == 1 | gamma == 0.8]
ggdat[, `:=`(n = factor(n),
             subsets = factor(subsets),
             method = car::recode(method, "'ps'='PS';'cbps'='CBPS';'svmLinear2'='SVM'"))]
ggdat[, `:=`(method = factor(method, levels = c('PS', 'CBPS', 'SVM')))]

ggplot(ggdat, aes(x = subsets, y = time_out)) +
  geom_boxplot() +
  facet_grid(method ~ n, scales = 'free_y', labeller = label_both) +
  xlab('Subsets') +
  ylab('Time elapsed (seconds)') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle = -25))
ggsave(file.path(image_path, 'time_comp.pdf'), width = 7, height = 8)

