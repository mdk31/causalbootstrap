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
  out <- cbind(y, Tr, X1, X2)
  out <- data.table::as.data.table(out)
  return(out)
}

# Non-independent
kangschafer3_dep <- function(n, te, beta_overlap = 0.5, sigma, rho = 0.8) {
  X <- mvtnorm::rmvnorm(n, sigma = diag(rho, nrow = 2, ncol = 2))
  X1 <- X[1, ]
  X2 <- X[2, ]
  prt <- 1/(1 + exp(-beta_overlap*(X1 + X2)))
  
  Tr <- rbinom(n, 1, prt)
  
  Y0  <- X1 + X2 + rnorm(n, 0, sigma)
  Y1  <- Y0 + te
  y   <- Y0*(1-Tr) + Y1*(Tr)
  out <- cbind(y, Tr, X1, X2)
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

make_weights_ML <- function(formula, data, method, normed = FALSE, ...){
  assertthat::assert_that(is.logical(normed))
  treat <- as.character(as.formula(formula)[2])
  if(!is.factor(data[[treat]])){
    data[[treat]] <- factor(data[[treat]])
  }
  vars <- all.vars(formula)
  vars <- vars[vars != treat]
  var_check <- sapply(vars, function(p) is.numeric(data[[p]]))
  var_check <- vars[var_check]
  if(length(var_check) > 0){
    for(i in var_check){
      data[[i]] <- scale(data[[i]])
    }
  }
  if(method == 'svm'){
    svm_out <- e1071::svm(formula = as.formula(formula),
                          data = data,
                          # kernel = 'linear',
                          type = 'C-classification',
                          # cost = 0.01,
                          probability = TRUE,
                          ...)
    preds <- attr(predict(svm_out, newdata = data, probability = TRUE),
                  "probabilities")[, "1"]
  } else if(method == 'ranger'){
    rang_out <- ranger::ranger(formula = as.formula(formula),
                               data = data,
                               probability = TRUE, ...)
    preds <- predict(rang_out, data = data)$predictions[, '1']
  }

  wts_svm <- 1/preds*(data[[treat]] == 1) + (1/(1-preds))*(data[[treat]] == 0)

  if(normed){
    ns <- abs(tapply(wts_svm, data[[treat]], sum))
    ns <- ns['1']*(data[[treat]] == 1) + (data[[treat]] == 0)*ns['0']
  } else{
    ns <- 1L
  }

  wts <- wts_svm/ns

  return(wts)
}

weight_blb <- function(data, R, pi_formula = NULL, Y, W, subsets,
                       type = 'obs', b = NULL, method = 'ps', ci_prep = 'perc', ...){
  stopifnot(Y %in% names(data))
  stopifnot(W %in% names(data))
  assertthat::assert_that(all(data[[W]] %in% 0:1))
  assertthat::assert_that(length(unique(data[[W]])) > 1)

  n <- nrow(data)
  n0 <- sum(data[[W]] == 0)
  n1 <- n - n0
  r <- round(R/subsets)
  # if(is.null(b)){
  #   obs_per_set <- n/subsets
  #   partition <- sample(seq(1, n, by = 1), n)
  # } else{
  #   assertthat::assert_that(is.integer(b))
  #   partition <- sample(seq(1, n, by = 1), b*subsets, replace = FALSE)
  #   obs_per_set <- b
  # }
  # partition <- split(partition, f = rep(1:subsets, each = obs_per_set))
  partition <- make_partition(n = n, subsets = subsets, b = b)
  lapply(partition, function(p){
    part_dat <- data[p]
    if(type == 'obs'){
      if(method %in% c('svm', 'ranger')){
        wts <- make_weights_ML(formula = as.formula(pi_formula),
                            data = part_dat,
                            method = method,
                            normed = TRUE, ...)
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
    ws <- split(wts, part_dat[[W]])
    ys <- split(part_dat[[Y]], part_dat[[W]])
    w0 <- ws[['0']]
    w1 <- ws[['1']]

    M1 <- rmultinom(r, n1, prob = w1)
    y1 <- colSums(ys[['1']]*M1)/n1
    M0 <- rmultinom(r, n0, prob = w0)
    y0 <- colSums(ys[['0']]*M0)/n0

    if(ci_prep == 'perc'){
      return(y1 - y0)
    } else{
      # Calculate subset theta and quantile sd
      theta <- sum(ys[['1']]*w1) - sum(ys[['0']]*w0)
      list(theta = theta, estim = y1 - y0)
    }
  })
}


