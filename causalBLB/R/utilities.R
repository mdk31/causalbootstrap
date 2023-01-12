
make_partition <- function(n, subsets, b = NULL){
  part_idx <- seq(1, n, by = 1)
  replicate(subsets, {
    sample(part_idx, size = b, replace = FALSE)
  }, simplify = FALSE)

}

# binary_ate <- function(W, Y, Y0 = NULL, Y1 = NULL){
#   chk1 <- (!is.null(W) & !is.null(Y))
#   chk2 <- (!is.null(Y0) & !is.null(Y1))
#   assertthat::assert_that(chk1 || chk2)
#   # if()
# }
