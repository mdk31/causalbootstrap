

causal_blb <- function(data, r, pi_formula = NULL, Y, W, subsets,
           type = 'obs', gamma = 0.8, method = 'ps', ci_type = 'perc', ...){

    stopifnot(Y %in% names(data))
    stopifnot(W %in% names(data))
    assertthat::assert_that(levels(factor(data[[W]]))[1] %in% c("FALSE", "0", 0) &
                              (levels(factor(data[[W]]))[2] %in% c("TRUE", "1", 1)) &
                              (length(levels(factor(data[[W]]))) == 2))
    treat <- factor(data[[W]])
    levels(treat) <- c(0, 1)

    n <- nrow(data)
    n0 <- sum(treat == 0)
    n1 <- n - n0
    b <- n^gamma
    partition <- make_partition(n = n, subsets = subsets, b = b)
    lapply(partition, function(p){
      part_dat <- data[p]
      if(type == 'obs'){
        wts <- make_weights()
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

}
