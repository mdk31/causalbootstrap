

make_weights <- function(formula, data, method = 'ps', ...){
  A <- list(...)

  if (is.null(formula) || is.null(class(formula))) {
    stop("'formula' must be a formula relating treatment to covariates.", call. = FALSE)
  }

  treat <- as.character(as.formula(formula)[2])
  A[['method']] <- method
  A[['formula']] <- formula
  A[['data']] <- data

  wo <- do.call('WeightIt::weightit', A)

  if(normed){
    ns <- abs(tapply(wo$weights, data[[treat]], sum))
    ns <- ns['1']*(data[[treat]] == 1) + (data[[treat]] == 0)*ns['0']
  } else{
    ns <- 1L
  }
  wts <- wo$weights/ns

  return(wts)

}
