

causal_blb <- function(data, r, pi_formula = NULL, Y, W, subsets, normed = TRUE,
           type = 'obs', gamma = 0.8, method = 'ps', ...){

    assertthat::assert_that(W %in% names(data), msg = 'Treatment variable W not in the data')
    assertthat::assert_that(Y %in% names(data), msg = 'Outcome variable Y not in the data')
    assertthat::assert_that(r > 0, msg = 'r must be greater than 0')
    treat <- factor(data[[W]])
    assertthat::assert_that(levels(treat)[1] %in% c("FALSE", "0", 0) &
                              (levels(treat)[2] %in% c("TRUE", "1", 1)) &
                              (length(levels(treat)) == 2))
    levels(treat) <- c(0, 1)
    n <- nrow(data)
    n0 <- sum(treat == 0)
    n1 <- n - n0
    b <- n^gamma

    partition <- make_partition(n = n, subsets = subsets, b = b)
    boot0 <- lapply(partition, function(p){
      part_dat <- data[p]
      if(type == 'obs'){
        wts <- make_weights(formula = as.formula(pi_formula),
                            data = part_dat,
                            method = method,
                            normed = normed)
      } else{
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

      theta <- sum(ys[['1']]*w1) - sum(ys[['0']]*w0)
      list(theta = theta, estim = y1 - y0)
    })

    thetas <- sapply(boot0, `[[`, 1)
    estims <- lapply(boot0, `[[`, 2)

    t <- do.call(rbind, estims)

    boot0 <- list(t = t,
                  t0 = thetas)

    class(boot0) <- 'causalblb'


}

print.causalblb <- function(x){

  cat('\nCAUSAL BOOTSTRAP\n\n')

}

boot.ci <- function(boot.out, conf = 0.95, type = 'all',  ...){

  call <- match.call()
  output <- list(call = call)

  if(any(type == 'all' | type == 'perc')){
    output <- c(output, list(percent = percblb.ci(boot.out$t, conf = conf)))
  }

  if(any(type == 'all' | type == 'norm')){
    output <- c(output, list(normal = normblb.ci(boot.out, conf = conf)))
  }

  class(output) <- "causalblbci"
  output

}

percblb.ci <- function(t, conf = 0.95){
  alpha <- (1 + c(-conf, conf))/2
  qq <- apply(t, 1, quantile, probs = alpha)
  qq <- Reduce(`+`, qq)/length(qq)
  cbind(conf, qq[1], qq[2])
}

normblb.ci <- function(boot.out, conf = 0.95){
  if (!is.null(boot.out)) t0 <- boot.out$t0
  else stop("bootstrap output object required")
  t <- boot.out$t

  merr <- apply(t, 1, sd)*qnorm((1+conf)/2)
  lower <- mean(t0 - merr)
  upper <- mean(t0 + merr)
  cbind(conf, lower, upper)
}
