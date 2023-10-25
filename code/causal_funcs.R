inv_logit <- function(x){
  1/(1 + exp(-x))
}


kangschafer3 <- function(n, te, beta_overlap = 0.5, sigma) {
  X1 <- rnorm(n, 0, 1)
  X2 <- rnorm(n, 0, 1)
  prt <- 1/(1 + exp(-beta_overlap*(X1 + X2)))

  Tr <- rbinom(n, 1, prt)

  Y0  <- X1 + X2 + rnorm(n, 0, sigma)
  Y1  <- Y0 + te
  y   <- Y0*(1-Tr) + Y1*(Tr)
  out <- cbind(y, Tr, X1, X2, Y1, Y0)
  out <- data.table::as.data.table(out)
  return(out)
}

# Non-independent
kangschafer3_dep <- function(n, te, beta_overlap = 0.5, sigma, rho = 0.8) {
  Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
  X <- mvtnorm::rmvnorm(n, sigma = Sigma)
  X1 <- X[, 1]
  X2 <- X[, 2]
  prt <- 1/(1 + exp(-beta_overlap*(X1 + X2)))
  
  Tr <- rbinom(n, 1, prt)
  
  Y0  <- X1 + X2 + rnorm(n, 0, sigma)
  Y1  <- Y0 + te
  y   <- Y0*(1-Tr) + Y1*(Tr)
  out <- cbind(y, Tr, X1, X2, Y1, Y0)
  out <- data.table::as.data.table(out)
  return(out)
}

# Multiple predictors
kangschafer3_mult <- function(n, te, beta_overlap = 0.1, sigma, preds = 10) {
  Sigma <- diag(preds)
  X <- mvtnorm::rmvnorm(n, sigma = Sigma)
  X_comb <- rowSums(X)
  prt <- 1/(1 + exp(-beta_overlap*(X_comb)))
  
  Tr <- rbinom(n, 1, prt)
  
  Y0  <- X_comb + rnorm(n, 0, sigma)
  Y1  <- Y0 + te
  y   <- Y0*(1-Tr) + Y1*(Tr)
  out <- cbind(y, Tr, X, Y1, Y0)
  out <- data.table::as.data.table(out)
  return(out)
}

# Heterogeneity
kangschafer3_het <- function(n, te, beta_overlap = 0.5, sigma) {
  X1 <- rnorm(n, 0, 1)
  X2 <- rnorm(n, 0, 1)
  prt <- 1/(1 + exp(-beta_overlap*(X1 + X2)))
  
  Tr <- rbinom(n, 1, prt)
  
  Y0  <- X1 + X2 + rnorm(n, 0, sigma)
  Y1  <- Y0 + te + X1
  y   <- Y0*(1-Tr) + Y1*(Tr)
  out <- cbind(y, Tr, X1, X2, Y0, Y1)
  out <- data.table::as.data.table(out)
  return(out)
}
# bias = \hat{ATE} - (Y1 - Y0)

# Misspecified logit model (Z1 and Z1 instead of X1 and X2)
kangschafer3_mis <- function(n, te, beta_overlap = 0.5, sigma) {
  X1 <- rnorm(n, 0, 1)
  X2 <- rnorm(n, 0, 1)
  
  Z1 <- X1^3/2
  Z2 <- 0.25*(X1 + X2)^2
  
  prt <- 1/(1 + exp(-beta_overlap*(Z1 + Z2)))
  
  Tr <- rbinom(n, 1, prt)
  
  Y0  <- Z1 + Z2 + rnorm(n, 0, sigma)
  Y1  <- Y0 + te + Z1
  y   <- Y0*(1-Tr) + Y1*(Tr)
  out <- cbind(y, Tr, X1, X2, Y0, Y1)
  out <- data.table::as.data.table(out)
  return(out)
}

make_partition <- function(n, subsets, b = NULL){
  part_idx <- seq(1, n, by = 1)
  replicate(subsets, {
      sample(part_idx, size = b, replace = FALSE)
    }, simplify = FALSE)

}

make_weights <- function(formula, data, method = 'ps', normed = FALSE){

  assertthat::assert_that(is.logical(normed))
  treat <- as.character(as.formula(formula)[2])
  wi <- WeightIt::weightit(formula = as.formula(formula),
                           data = data,
                           method = method,
                           estimand = 'ATE',
                           over = FALSE)
  if(normed){
    ns <- abs(tapply(wi$weights, data[[treat]], sum))
    ns <- ns['1']*(data[[treat]] == 1) + (1 - data[[treat]] == 1)*ns['0']
  } else{
    ns <- 1L
  }
  wts <- wi$weights/ns

  return(wts)
}

make_weights_ML <- function(formula, data, method, normed = FALSE, prop_grid, ...){
  assertthat::assert_that(is.logical(normed))
  treat <- as.character(as.formula(formula)[2])
  # Change to factor variable and change names to make compatible with caret
  data$Tr2 <- factor(data[[treat]])
  levels(data$Tr2) <- make.names(levels(data$Tr2))
  # Update prop_formula
  pi_formula <- update.formula(as.formula(formula), Tr2 ~ .)
  
  vars <- all.vars(formula)
  vars <- vars[vars != treat]
  var_check <- sapply(vars, function(p) is.numeric(data[[p]]))
  var_check <- vars[var_check]
  if(length(var_check) > 0){
    for(i in var_check){
      data[[i]] <- scale(data[[i]])
    }
  }
  prop_mod <- caret::train(form = pi_formula,
                           data = data,
                           method = method,
                           type = 'C-classification',
                           probability = TRUE,
                           trControl = caret::trainControl(method = 'none', classProbs = TRUE),
                           tuneGrid = prop_grid,
                           verbose = FALSE,
                           scale = FALSE,
                           ...)
  ps <- predict(prop_mod, data, 'prob')[, 2]
  wts <- 1/ps*(data[[treat]] == 1) + (1/(1-ps))*(data[[treat]] == 0)
  if(normed){
    ns <- abs(tapply(wts, data[[treat]], sum))
    ns <- ns['1']*(data[[treat]] == 1) + (data[[treat]] == 0)*ns['0']
  } else{
    ns <- 1L
  }
  wts <- wts/ns
  return(wts)
}

# make_weights_ML <- function(formula, data, method, normed = FALSE, ...){
#   assertthat::assert_that(is.logical(normed))
#   treat <- as.character(as.formula(formula)[2])
#   if(!is.factor(data[[treat]])){
#     data[[treat]] <- factor(data[[treat]])
#   }
#   vars <- all.vars(formula)
#   vars <- vars[vars != treat]
#   var_check <- sapply(vars, function(p) is.numeric(data[[p]]))
#   var_check <- vars[var_check]
#   if(length(var_check) > 0){
#     for(i in var_check){
#       data[[i]] <- scale(data[[i]])
#     }
#   }
#   if(method == 'svm'){
#     svm_out <- e1071::svm(formula = as.formula(formula),
#                           data = data,
#                           type = 'C-classification',
#                           probability = TRUE,
#                           ...)
#     preds <- attr(predict(svm_out, newdata = data, probability = TRUE),
#                   "probabilities")[, "1"]
#   } else if(method == 'ranger'){
#     rang_out <- ranger::ranger(formula = as.formula(formula),
#                                data = data,
#                                probability = TRUE, ...)
#     preds <- predict(rang_out, data = data)$predictions[, '1']
#   }
# 
#   wts_svm <- 1/preds*(data[[treat]] == 1) + (1/(1-preds))*(data[[treat]] == 0)
# 
#   if(normed){
#     ns <- abs(tapply(wts_svm, data[[treat]], sum))
#     ns <- ns['1']*(data[[treat]] == 1) + (data[[treat]] == 0)*ns['0']
#   } else{
#     ns <- 1L
#   }
# 
#   wts <- wts_svm/ns
# 
#   return(wts)
# }

weight_blb <- function(data, R, pi_formula = NULL, Y, W, subsets,
                       type = 'obs', b = NULL, method = 'ps', ci_prep = 'perc', prop_grid = NULL, ...){
  stopifnot(Y %in% names(data))
  stopifnot(W %in% names(data))
  assertthat::assert_that(all(data[[W]] %in% 0:1))
  assertthat::assert_that(length(unique(data[[W]])) > 1)

  n <- nrow(data)
  n0 <- sum(data[[W]] == 0)
  n1 <- n - n0

  true_ATE <- mean(data$Y1) - mean(data$Y0)
  partition <- make_partition(n = n, subsets = subsets, b = b)
  lapply(partition, function(p){
    part_dat <- data[p]
    if(type == 'obs'){
      if(!method %in% c('ps', 'cbps')){
        wts <- make_weights_ML(formula = as.formula(pi_formula),
                            data = part_dat,
                            method = method,
                            normed = TRUE,
                            prop_grid = prop_grid, 
                            ...)
      } else{
        wts <- make_weights(formula = as.formula(pi_formula),
                            data = part_dat,
                            method = method,
                            normed = TRUE)
      }
    } else{
      # wts <- with(part_dat, W/n1 + (1 - W)/n0)*(1/2)
      # Normed so weights sum to 1 in each group
      wts <- (part_dat[[W]]/n1 + (1- part_dat[[W]])/n0)
    }
    theta0 <- sum(part_dat[[W]]*part_dat[[Y]]*wts - (1-part_dat[[W]])*part_dat[[Y]]*wts)
    ws <- split(wts, part_dat[[W]])
    ys <- split(part_dat[[Y]], part_dat[[W]])
    w0 <- ws[['0']]
    w1 <- ws[['1']]

    M1 <- rmultinom(R, n1, prob = w1)
    y1 <- colSums(ys[['1']]*M1)/n1
    M0 <- rmultinom(R, n0, prob = w0)
    y0 <- colSums(ys[['0']]*M0)/n0

    if(ci_prep == 'perc'){
      lst <- list(theta0 = theta0, theta_reps = y1 - y0, truth = true_ATE)
      return(lst)
    } else{
      # Calculate subset theta and quantile sd
      theta <- sum(ys[['1']]*w1) - sum(ys[['0']]*w0)
      list(theta = theta, estim = y1 - y0)
    }
  })
}


