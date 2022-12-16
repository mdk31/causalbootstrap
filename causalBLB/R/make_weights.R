

make_weights <- function(formula, method = 'ps'){
  A <- list(...)

  if (is.null(formula) || is.null(class(formula))) {
    stop("'formula' must be a formula relating treatment to covariates.", call. = FALSE)
  }

  treat <- as.character(as.formula(formula)[2])

  wo <- do.call('WeightIt::weightit', A)

  wi <- WeightIt::weightit(formula = formula,
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
